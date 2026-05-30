import sqlite3
from pathlib import Path
from flask import Blueprint, render_template, jsonify, request, g, redirect
import math

timecourse_bp = Blueprint("timecourse", __name__, url_prefix="/timecourse")

DB_PATH = str(Path(__file__).resolve().parent / "timecourse.db")

# Clusters excluded from all displays (temporary exclusion list)
HIDDEN_GENE_CLUSTERS = frozenset({'gene_cluster_7'})
HIDDEN_PEAK_CLUSTERS = frozenset({'peak_cluster_5'})
_HIDDEN_NODES = HIDDEN_GENE_CLUSTERS | HIDDEN_PEAK_CLUSTERS


def get_db():
    if "tc_db" not in g:
        g.tc_db = sqlite3.connect(DB_PATH)
        g.tc_db.row_factory = sqlite3.Row
    return g.tc_db


@timecourse_bp.teardown_request
def close_db(exc):
    db = g.pop("tc_db", None)
    if db is not None:
        db.close()


def rows_to_dicts(rows):
    return [dict(r) for r in rows]


# ── Page routes ───────────────────────────────────────────────────────────────

def _gc_list(db):
    return [r[0] for r in db.execute(
        "SELECT DISTINCT cluster FROM gene_cluster_profiles ORDER BY cluster")
        if r[0] not in HIDDEN_GENE_CLUSTERS]

def _pc_list(db):
    return [r[0] for r in db.execute(
        "SELECT DISTINCT cluster FROM peak_cluster_profiles ORDER BY cluster")
        if r[0] not in HIDDEN_PEAK_CLUSTERS]


def _chip_paths(profile, w, h):
    """Return (line_path, area_path) SVG path strings for a sidebar chip sparkline."""
    vals = [p["val"] for p in profile]
    n = len(vals)
    if n < 2:
        return "", ""
    xs = [3 + (i / (n - 1)) * (w - 6) for i in range(n)]
    mn, mx = min(vals), max(vals)
    rng = mx - mn if mx != mn else 1.0
    pad = 2
    ys = [h - pad - ((v - mn) / rng) * (h - 2 * pad) for v in vals]
    pts = " L ".join(f"{x:.1f},{y:.1f}" for x, y in zip(xs, ys))
    return f"M {pts}", f"M {pts} L {xs[-1]:.1f},{h} L {xs[0]:.1f},{h} Z"


def _chip_data(db, table, hidden, gc=True):
    """Return cluster list sorted by temporal peak, with HSL gradient colors and sparkline paths."""
    profiles = {}
    for r in db.execute(
            f"SELECT cluster, timepoint, value FROM {table} ORDER BY cluster, rowid"):
        if r[0] not in hidden:
            profiles.setdefault(r[0], []).append({"tp": r[1], "val": r[2]})
    result = []
    for cid, prof in profiles.items():
        vals = [p["val"] for p in prof]
        peak_idx = vals.index(max(vals))
        line_path, area_path = _chip_paths(prof, 74, 22)
        result.append({
            "id": cid,
            "label": cid.replace("gene_cluster_", "GC").replace("peak_cluster_", "PC"),
            "peak_idx": peak_idx,
            "peak_tp": prof[peak_idx]["tp"],
            "line_path": line_path,
            "area_path": area_path,
        })
    result.sort(key=lambda c: c["peak_idx"])
    n = len(result)
    for i, c in enumerate(result):
        frac = i / max(1, n - 1)
        if gc:
            L = round(82 - frac * 47)
            c["bg_color"] = f"hsl(18,68%,{L}%)"
            c["label_color"] = f"hsl(18,68%,{min(L, 44)}%)"
        else:
            L = round(85 - frac * 57)
            c["bg_color"] = f"hsl(214,58%,{L}%)"
            c["label_color"] = f"hsl(214,58%,{min(L, 38)}%)"
    return result


@timecourse_bp.route("/")
def index():
    db = get_db()
    return render_template("timecourse/index.html",
                           gene_clusters=_gc_list(db),
                           peak_clusters=_pc_list(db),
                           gene_cluster_data=_chip_data(db, "gene_cluster_profiles",
                                                        HIDDEN_GENE_CLUSTERS, gc=True),
                           peak_cluster_data=_chip_data(db, "peak_cluster_profiles",
                                                        HIDDEN_PEAK_CLUSTERS, gc=False))


@timecourse_bp.route("/gene-cluster/<cluster_id>")
def gene_cluster(cluster_id):
    n = cluster_id.replace('gene_cluster_', '')
    return redirect(f'/perturbseq/module/GC{n}', 301)


@timecourse_bp.route("/peak-cluster/<cluster_id>")
def peak_cluster(cluster_id):
    n = cluster_id.replace('peak_cluster_', '')
    return redirect(f'/perturbseq/module/PC{n}', 301)


@timecourse_bp.route("/toggle")
def toggle():
    db = get_db()
    return render_template("timecourse/toggle.html",
                           gene_clusters=_gc_list(db),
                           peak_clusters=_pc_list(db))


@timecourse_bp.route("/dual")
def dual():
    db = get_db()
    return render_template("timecourse/dual.html",
                           gene_clusters=_gc_list(db),
                           peak_clusters=_pc_list(db))


@timecourse_bp.route("/clusters")
def clusters():
    db = get_db()

    def build_cluster_list(table, id_col, color_map=None, default_color="#e07b39"):
        rows = db.execute(
            f"SELECT cluster, timepoint, value FROM {table} ORDER BY cluster, rowid"
        ).fetchall()
        by_cluster = {}
        for r in rows:
            cid = r[0]
            if cid not in by_cluster:
                by_cluster[cid] = []
            by_cluster[cid].append({"tp": r[1], "val": r[2]})

        hidden = HIDDEN_GENE_CLUSTERS if table == "gene_cluster_profiles" else HIDDEN_PEAK_CLUSTERS
        result = []
        for cid, profile in sorted(by_cluster.items()):
            if cid in hidden:
                continue
            vals = [p["val"] for p in profile]
            peak_idx = vals.index(max(vals))
            peak_tp = profile[peak_idx]["tp"]
            color = (color_map.get(cid, "#aaa") if color_map else default_color)
            svg_path = _sparkline_path(vals, 150, 55)
            result.append({
                "id": cid,
                "profile": profile,
                "peak_tp": peak_tp,
                "color": color,
                "svg_path": svg_path,
            })
        result.sort(key=lambda c: [p["val"] for p in c["profile"]].index(
            max(p["val"] for p in c["profile"])
        ))
        return result

    gc_colors = {
        "gene_cluster_1": "#4E79A7", "gene_cluster_2": "#F28E2B",
        "gene_cluster_3": "#59A14F", "gene_cluster_4": "#E15759",
        "gene_cluster_5": "#B6992D", "gene_cluster_6": "#499894",
        "gene_cluster_7": "#79706E",
    }

    gene_clusters_data = build_cluster_list("gene_cluster_profiles", "cluster", gc_colors)
    peak_clusters_data = build_cluster_list("peak_cluster_profiles", "cluster")

    return render_template("timecourse/clusters.html",
                           gene_clusters=_gc_list(db),
                           peak_clusters=_pc_list(db),
                           gene_clusters_data=gene_clusters_data,
                           peak_clusters_data=peak_clusters_data)


def _sparkline_path(values, width, height):
    if len(values) < 2:
        return ""
    xs = [(i / (len(values) - 1)) * width for i in range(len(values))]
    mn, mx = min(values), max(values)
    rng = mx - mn if mx != mn else 1.0
    ys = [height - ((v - mn) / rng) * height for v in values]
    pts = " L ".join(f"{x:.1f},{y:.1f}" for x, y in zip(xs, ys))
    return f"M {pts}"


# ── API routes ────────────────────────────────────────────────────────────────

@timecourse_bp.route("/api/network")
def api_network():
    db = get_db()
    nes_padj      = float(request.args.get("nes_padj",      0.05))
    binding_fdr   = float(request.args.get("binding_fdr",   0.001))
    min_or        = float(request.args.get("min_or",         5.0))
    min_pg_or     = float(request.args.get("min_pg_or",      1.0))
    binding_src   = request.args.get("binding_source", "footprint")  # chipseq|footprint|both
    show_pert     = request.args.get("show_perturbation", "true").lower() != "false"
    show_interact = request.args.get("show_interaction",   "true").lower() != "false"

    # Pre-fetch profiles + compute gradient rank fractions for cluster nodes
    gc_data, pc_data = {}, {}
    for r in db.execute(
            "SELECT cluster, timepoint, value FROM gene_cluster_profiles ORDER BY cluster, rowid"):
        if r[0] not in HIDDEN_GENE_CLUSTERS:
            entry = gc_data.setdefault(r[0], {"vals": [], "profile": []})
            entry["vals"].append(r[2])
            entry["profile"].append({"timepoint": r[1], "value": round(r[2], 4)})
    for r in db.execute(
            "SELECT cluster, timepoint, value FROM peak_cluster_profiles ORDER BY cluster, rowid"):
        if r[0] not in HIDDEN_PEAK_CLUSTERS:
            entry = pc_data.setdefault(r[0], {"vals": [], "profile": []})
            entry["vals"].append(r[2])
            entry["profile"].append({"timepoint": r[1], "value": round(r[2], 4)})

    def _pidx(entry): return entry["vals"].index(max(entry["vals"])) if entry["vals"] else 0
    gc_sorted = sorted(gc_data.keys(), key=lambda c: _pidx(gc_data[c]))
    pc_sorted = sorted(pc_data.keys(), key=lambda c: _pidx(pc_data[c]))
    gc_rank_frac = {c: round(i / max(1, len(gc_sorted) - 1), 4)
                   for i, c in enumerate(gc_sorted)}
    pc_rank_frac = {c: round(i / max(1, len(pc_sorted) - 1), 4)
                   for i, c in enumerate(pc_sorted)}

    nodes = {}
    edges = []

    def add_node(node_id, node_type):
        if node_id not in nodes:
            node_data = {"id": node_id, "type": node_type}
            if node_type == "gene_cluster" and node_id in gc_data:
                node_data["peak_rank_frac"] = gc_rank_frac.get(node_id, 0)
                node_data["profile"] = gc_data[node_id]["profile"]
            elif node_type == "peak_cluster" and node_id in pc_data:
                node_data["peak_rank_frac"] = pc_rank_frac.get(node_id, 0)
                node_data["profile"] = pc_data[node_id]["profile"]
            nodes[node_id] = node_data

    # Gene cluster and peak cluster nodes always present (excluding hidden)
    for r in db.execute("SELECT DISTINCT cluster FROM gene_cluster_profiles"):
        if r[0] not in HIDDEN_GENE_CLUSTERS:
            add_node(r[0], "gene_cluster")
    for r in db.execute("SELECT DISTINCT cluster FROM peak_cluster_profiles"):
        if r[0] not in HIDDEN_PEAK_CLUSTERS:
            add_node(r[0], "peak_cluster")

    # Perturbation layer: TF → gene_cluster
    if show_pert:
        for r in db.execute(
            "SELECT tf, gene_cluster, nes, padj FROM perturbation_edges WHERE padj <= ?",
            (nes_padj,)):
            if r["gene_cluster"] in HIDDEN_GENE_CLUSTERS:
                continue
            add_node(r["tf"], "tf")
            edges.append({
                "source": r["tf"], "target": r["gene_cluster"],
                "type": "perturbation", "nes": round(r["nes"], 3), "padj": r["padj"],
            })

    # Interaction layer: TF → peak_cluster and peak_cluster → gene_cluster
    if show_interact:
        src_filter = "" if binding_src == "both" else f" AND source='{binding_src}'"
        for r in db.execute(
            f"SELECT tf, peak_cluster, odds_ratio, fdr, source "
            f"FROM binding_edges WHERE fdr <= ? AND odds_ratio >= ?{src_filter}",
            (binding_fdr, min_or)):
            if r["peak_cluster"] in HIDDEN_PEAK_CLUSTERS:
                continue
            add_node(r["tf"], "tf")
            edges.append({
                "source": r["tf"], "target": r["peak_cluster"],
                "type": "binding", "source_type": r["source"],
                "odds_ratio": round(r["odds_ratio"], 3), "fdr": r["fdr"],
            })
        # Peak → gene edges (filtered by min odds ratio)
        for r in db.execute(
            "SELECT peak_cluster, gene_cluster, odds_ratio, fdr FROM peak_gene_edges "
            "WHERE odds_ratio >= ?", (min_pg_or,)):
            if r["peak_cluster"] in HIDDEN_PEAK_CLUSTERS or r["gene_cluster"] in HIDDEN_GENE_CLUSTERS:
                continue
            edges.append({
                "source": r["peak_cluster"], "target": r["gene_cluster"],
                "type": "peak_gene",
                "odds_ratio": round(r["odds_ratio"], 3), "fdr": r["fdr"],
            })

    return jsonify({"nodes": list(nodes.values()), "edges": edges})


@timecourse_bp.route("/api/cluster/gene/<cluster_id>")
def api_gene_cluster(cluster_id):
    db = get_db()
    profile = rows_to_dicts(db.execute(
        "SELECT timepoint, value FROM gene_cluster_profiles WHERE cluster=? ORDER BY rowid",
        (cluster_id,)).fetchall())
    top_tfs = rows_to_dicts(db.execute(
        "SELECT tf, nes, padj FROM perturbation_edges WHERE gene_cluster=? "
        "ORDER BY abs(nes) DESC LIMIT 30",
        (cluster_id,)).fetchall())
    assoc_peaks = rows_to_dicts(db.execute(
        "SELECT peak_cluster, odds_ratio, fdr FROM peak_gene_edges "
        "WHERE gene_cluster=? ORDER BY odds_ratio DESC",
        (cluster_id,)).fetchall())
    return jsonify({"profile": profile, "top_tfs": top_tfs, "assoc_peaks": assoc_peaks})


@timecourse_bp.route("/api/cluster/peak/<cluster_id>")
def api_peak_cluster(cluster_id):
    db = get_db()
    profile = rows_to_dicts(db.execute(
        "SELECT timepoint, value FROM peak_cluster_profiles WHERE cluster=? ORDER BY rowid",
        (cluster_id,)).fetchall())
    top_tfs = rows_to_dicts(db.execute(
        "SELECT tf, odds_ratio, fdr, source FROM binding_edges "
        "WHERE peak_cluster=? ORDER BY odds_ratio DESC LIMIT 30",
        (cluster_id,)).fetchall())
    assoc_genes = rows_to_dicts(db.execute(
        "SELECT gene_cluster, odds_ratio, fdr FROM peak_gene_edges "
        "WHERE peak_cluster=? ORDER BY odds_ratio DESC",
        (cluster_id,)).fetchall())
    return jsonify({"profile": profile, "top_tfs": top_tfs, "assoc_genes": assoc_genes})


@timecourse_bp.route("/sankey")
def sankey():
    db = get_db()
    return render_template("timecourse/sankey.html",
                           gene_clusters=_gc_list(db),
                           peak_clusters=_pc_list(db))


@timecourse_bp.route("/api/sankey")
def api_sankey():
    db = get_db()
    binding_fdr = float(request.args.get("binding_fdr", 0.001))
    min_or      = float(request.args.get("min_or", 5.0))
    nes_padj    = float(request.args.get("nes_padj", 0.05))
    binding_src = request.args.get("binding_source", "footprint")

    nodes = []
    node_index = {}

    def add_node(node_id, node_type):
        if node_id not in node_index:
            node_index[node_id] = len(nodes)
            nodes.append({"id": node_id, "type": node_type})
        return node_index[node_id]

    # Cluster nodes always present (excluding hidden)
    for r in db.execute("SELECT DISTINCT cluster FROM peak_cluster_profiles ORDER BY cluster"):
        if r[0] not in HIDDEN_PEAK_CLUSTERS:
            add_node(r[0], "peak_cluster")
    for r in db.execute("SELECT DISTINCT cluster FROM gene_cluster_profiles ORDER BY cluster"):
        if r[0] not in HIDDEN_GENE_CLUSTERS:
            add_node(r[0], "gene_cluster")

    links = []

    src_filter = "" if binding_src == "both" else f" AND source='{binding_src}'"
    for r in db.execute(
        f"SELECT tf, peak_cluster, odds_ratio FROM binding_edges "
        f"WHERE fdr <= ? AND odds_ratio >= ?{src_filter}", (binding_fdr, min_or)):
        if r[1] in HIDDEN_PEAK_CLUSTERS:
            continue
        s = add_node(r[0], "tf")
        t = node_index[r[1]]
        links.append({"source": s, "target": t, "value": r[2], "layer": "binding"})

    for r in db.execute("SELECT peak_cluster, gene_cluster, odds_ratio FROM peak_gene_edges"):
        if r[0] in HIDDEN_PEAK_CLUSTERS or r[1] in HIDDEN_GENE_CLUSTERS:
            continue
        links.append({
            "source": node_index[r[0]],
            "target": node_index[r[1]],
            "value": r[2],
            "layer": "peak_gene",
        })

    for r in db.execute(
        "SELECT tf, gene_cluster, nes, padj FROM perturbation_edges WHERE padj <= ?",
        (nes_padj,)):
        if r[1] in HIDDEN_GENE_CLUSTERS:
            continue
        s = add_node(r[0], "tf")
        t = node_index[r[1]]
        links.append({
            "source": s, "target": t,
            "value": abs(r[2]), "nes": r[2],
            "layer": "perturbation",
        })

    return jsonify({"nodes": nodes, "links": links})
