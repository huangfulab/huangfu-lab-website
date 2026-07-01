import csv
import json
import math
import re
import sqlite3
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.stats import gaussian_kde
from flask import Blueprint, render_template, jsonify, request, g, redirect, url_for, abort, send_from_directory

PERTURBSEQ_PREFIX = '/endoderm-perturbseq'

perturbseq_bp = Blueprint('perturbseq', __name__, url_prefix=PERTURBSEQ_PREFIX)

DB_PATH      = str(Path(__file__).resolve().parent.parent / "data" / "db" / "tf-perturbseq-v6.db")

def _last_commit_ts() -> int | None:
    try:
        out = subprocess.check_output(
            ['git', 'log', '-1', '--format=%ct'],
            cwd=str(Path(__file__).resolve().parent.parent),
            stderr=subprocess.DEVNULL, text=True,
        ).strip()
        return int(out) if out else None
    except Exception:
        return None

def _fmt_time_ago(ts: int | None) -> str | None:
    if ts is None:
        return None
    diff = int(datetime.now(timezone.utc).timestamp()) - ts
    if diff < 60:
        n = diff
        return f"{n} second{'s' if n != 1 else ''} ago"
    if diff < 3600:
        n = diff // 60
        return f"{n} minute{'s' if n != 1 else ''} ago"
    if diff < 86400:
        n = diff // 3600
        return f"{n} hour{'s' if n != 1 else ''} ago"
    if diff < 172800:
        return "yesterday"
    return datetime.fromtimestamp(ts, tz=timezone.utc).strftime("%-d %b %Y")

_LAST_COMMIT_TS: int | None = _last_commit_ts()
BW_DIR       = str(Path(__file__).resolve().parent.parent / "data" / "bw" / "atac")
BW_EXTRA_DIR = BW_DIR
BW_RNA_DIR   = str(Path(__file__).resolve().parent.parent / "data" / "bw" / "rna")
GTF_DIR      = str(Path(__file__).resolve().parent.parent / "data" / "gtf")

TC_CLUSTER_COLORS = {
    "gene_cluster_1": "#4E79A7",
    "gene_cluster_2": "#F28E2B",
    "gene_cluster_3": "#59A14F",
    "gene_cluster_4": "#E15759",
    "gene_cluster_5": "#B6992D",
    "gene_cluster_6": "#499894",
}
HIDDEN_TC = frozenset({'gene_cluster_7', 'peak_cluster_5'})
MODULE_COLORS = {
    'DE-1':'#4E79A7','DE-2':'#A0CBE8','DE-3':'#F28E2B',
    'DE-4':'#FFBE7D','DE-5':'#59A14F','DE-6':'#8CD17D','DE-7':'#B6992D',
    'DE-8':'#F1CE63','DE-9':'#499894','DE-10':'#86BCB6','DE-11':'#E15759',
    'DE-12':'#FF9D9A',
}
_TP_ORDER = ['ES_0h', 'DE_12h', 'DE_24h', 'DE_36h', 'DE_48h', 'DE_60h', 'DE_72h']
_SAMPLE_TIMEPOINT_RE = re.compile(r'_\d+$')

def _top_genes_by_pub(genes: list, n: int = 8) -> list:
    pub_counts = _load_gene_pub_counts()
    return sorted(genes, key=lambda g: pub_counts.get(g, 0), reverse=True)[:n]


def _genes_to_tpm_by_tp(db, genes: list) -> dict:
    """Map gene names → {timepoint → mean_tpm} by joining gene_table + bulk_expression.

    Strips replicate suffix from sample names (DE_12h_1 → DE_12h) to aggregate.
    """
    if not genes:
        return {}
    ph = ','.join('?' * len(genes))
    id_rows = db.execute(
        f'SELECT gene_id, gene_name FROM gene_table WHERE gene_name IN ({ph})', genes
    ).fetchall()
    gene_id_map = {r['gene_name']: r['gene_id'] for r in id_rows}
    gene_ids = list(gene_id_map.values())
    if not gene_ids:
        return {}
    ph2 = ','.join('?' * len(gene_ids))
    expr_rows = db.execute(
        f'SELECT gene_id, sample_name AS sample, tpm FROM bulk_expression WHERE gene_id IN ({ph2})', gene_ids
    ).fetchall()
    id_to_name = {v: k for k, v in gene_id_map.items()}
    raw: dict = defaultdict(lambda: defaultdict(list))
    for gene_id, sample, tpm in expr_rows:
        if tpm is not None:
            raw[gene_id][_SAMPLE_TIMEPOINT_RE.sub('', sample)].append(tpm)
    result: dict = {}
    for gid, by_tp in raw.items():
        result[id_to_name[gid]] = {tp: sum(v) / len(v) for tp, v in by_tp.items()}
    return result


def _compute_mean_zscore_profile(db, genes: list, timepoints: list) -> list:
    """Z-score each gene's trajectory (row-normalise), return mean ± SEM per timepoint."""
    if not genes or not timepoints:
        return []
    gene_tpm = _genes_to_tpm_by_tp(db, genes)
    tp_zscores: dict = defaultdict(list)
    for gene, by_tp in gene_tpm.items():
        vals = [by_tp.get(tp) for tp in timepoints]
        if None in vals:
            continue
        mu = sum(vals) / len(vals)
        sigma = math.sqrt(sum((v - mu) ** 2 for v in vals) / max(len(vals) - 1, 1))
        if sigma < 1e-6:
            continue
        for tp, v in zip(timepoints, vals):
            tp_zscores[tp].append((v - mu) / sigma)
    result = []
    for tp in timepoints:
        zs = tp_zscores.get(tp, [])
        if not zs:
            result.append({'timepoint': tp, 'mean_z': 0.0, 'sem_z': 0.0, 'sd_z': 0.0})
            continue
        mu_z = sum(zs) / len(zs)
        sd_z = math.sqrt(sum((z - mu_z) ** 2 for z in zs) / max(len(zs) - 1, 1))
        sem_z = sd_z / math.sqrt(len(zs))
        result.append({'timepoint': tp, 'mean_z': round(mu_z, 4), 'sem_z': round(sem_z, 4), 'sd_z': round(sd_z, 4)})
    return result


def _compute_mean_tpm_profile(db, genes: list, timepoints: list) -> list:
    """Return mean ± SD TPM per timepoint across genes."""
    if not genes or not timepoints:
        return []
    gene_tpm = _genes_to_tpm_by_tp(db, genes)
    tp_vals: dict = defaultdict(list)
    for by_tp in gene_tpm.values():
        for tp, v in by_tp.items():
            if tp in timepoints:
                tp_vals[tp].append(v)
    result = []
    for tp in timepoints:
        vals = tp_vals.get(tp, [])
        n = len(vals)
        if n == 0:
            result.append({'timepoint': tp, 'mean_tpm': 0.0, 'sd_tpm': 0.0, 'n_genes': 0})
        else:
            mu = sum(vals) / n
            sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / max(n - 1, 1))
            result.append({'timepoint': tp, 'mean_tpm': round(mu, 4), 'sd_tpm': round(sd, 4), 'n_genes': n})
    return result


# ── New-schema helpers (TfDatabase) ──────────────────────────────────────────

def resolve_gene(db, name: str) -> dict | None:
    """Resolve a display gene name to a gene_table row, falling back to gene_synonym.

    Returns dict with gene_id, gene_name, chr, start, end or None if not found.
    Used as the entry point for all new-schema gene lookups.
    """
    row = db.execute(
        'SELECT gene_id, gene_name, chr, chrom_start AS start, chrom_end AS end FROM gene_table WHERE gene_name=? LIMIT 1',
        (name,)
    ).fetchone()
    if not row:
        row = db.execute(
            'SELECT g.gene_id, g.gene_name, g.chr, g.chrom_start AS start, g.chrom_end AS end '
            'FROM gene_synonym s JOIN gene_table g USING(gene_id) WHERE s.synonym=? LIMIT 1',
            (name,)
        ).fetchone()
    return dict(row) if row else None


def _gene_expr_by_timepoint(db, gene_id: str) -> list:
    """Aggregate bulk_expression per timepoint for a single gene (new schema).

    Sample names like 'DE_12h_2' → timepoint 'DE_12h' by stripping trailing _N.
    Returns [{timepoint, mean_tpm, replicates}] in _TP_ORDER order.
    """
    rows = db.execute(
        "SELECT sample_name AS sample, tpm FROM bulk_expression WHERE gene_id=?", (gene_id,)
    ).fetchall()
    by_tp: dict = defaultdict(list)
    for sample, tpm in rows:
        if tpm is not None:
            tp = _SAMPLE_TIMEPOINT_RE.sub('', sample)
            by_tp[tp].append(tpm)
    result = []
    for tp in _TP_ORDER:
        vals = by_tp.get(tp, [])
        mean_tpm = round(sum(vals) / len(vals), 4) if vals else 0.0
        result.append({
            'timepoint': tp,
            'mean_tpm': mean_tpm,
            'replicates': [round(v, 4) for v in vals],
        })
    return result


def _module_expr_by_timepoint(db, module_id: int) -> list:
    """Compute mean ± SD TPM per timepoint across all genes in a module (new schema).

    Fetches gene_ids from gene_module_table, then aggregates bulk_expression.
    Returns [{timepoint, mean_tpm, sd_tpm, n_genes}] in _TP_ORDER order.
    """
    gene_ids = [r[0] for r in db.execute(
        "SELECT gene_id FROM gene_module_table WHERE module_id=? AND gene_id IS NOT NULL",
        (module_id,)
    ).fetchall()]
    if not gene_ids:
        return [{'timepoint': tp, 'mean_tpm': 0.0, 'sd_tpm': 0.0, 'n_genes': 0} for tp in _TP_ORDER]
    ph = ','.join('?' * len(gene_ids))
    rows = db.execute(
        f"SELECT sample_name AS sample, tpm FROM bulk_expression WHERE gene_id IN ({ph})", gene_ids
    ).fetchall()
    by_tp: dict = defaultdict(list)
    for sample, tpm in rows:
        if tpm is not None:
            tp = _SAMPLE_TIMEPOINT_RE.sub('', sample)
            by_tp[tp].append(tpm)
    result = []
    for tp in _TP_ORDER:
        vals = by_tp.get(tp, [])
        n = len(vals)
        if n == 0:
            result.append({'timepoint': tp, 'mean_tpm': 0.0, 'sd_tpm': 0.0, 'n_genes': 0})
        else:
            mu = sum(vals) / n
            sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / max(n - 1, 1))
            result.append({'timepoint': tp, 'mean_tpm': round(mu, 4), 'sd_tpm': round(sd, 4), 'n_genes': n})
    return result


# ─────────────────────────────────────────────────────────────────────────────

def _tc_display(cid: str) -> str:
    if cid.startswith('gene_cluster_'):
        return 'GC' + cid[len('gene_cluster_'):]
    if cid.startswith('peak_cluster_'):
        return 'PC' + cid[len('peak_cluster_'):]
    return cid

def _mfuzz_to_internal(name: str) -> str:
    """Convert 'cluster_N' (new DB module_name) to 'gene_cluster_N' (legacy internal ID)."""
    return 'gene_cluster_' + name.split('_')[1]

_expressed_tfs: frozenset | None = None
_lambert_tfs: frozenset | None = None
_perturbed_tfs: frozenset | None = None
_gene_pub_counts: dict | None = None

def _query_db_raw(sql: str, params=()) -> list:
    conn = sqlite3.connect(DB_PATH)
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()


def _load_lambert_tfs() -> frozenset:
    global _lambert_tfs
    if _lambert_tfs is None:
        try:
            path = Path(__file__).resolve().parent.parent.parent.parent / "data" / "tables" / "lambert-TF_table.csv"
            with open(path, newline='') as f:
                _lambert_tfs = frozenset(
                    row['Name'] for row in csv.DictReader(f)
                    if row.get('is_TF') == 'Yes'
                )
        except Exception:
            _lambert_tfs = frozenset()
    return _lambert_tfs

def _load_gene_pub_counts() -> dict:
    global _gene_pub_counts
    if _gene_pub_counts is None:
        try:
            rows = _query_db_raw("""
                SELECT gt.gene_name, pc.publication_count
                FROM gene_publication_count pc
                JOIN gene_table gt ON pc.gene_id = gt.gene_id
            """)
            _gene_pub_counts = {r[0]: r[1] for r in rows}
        except Exception:
            _gene_pub_counts = {}
    return _gene_pub_counts


def _load_expressed_tfs() -> frozenset:
    """TFs present in the perturbation library (in_perturbation_library=1 in gene_table)."""
    global _expressed_tfs
    if _expressed_tfs is None:
        try:
            rows = _query_db_raw(
                "SELECT gene_name FROM gene_table WHERE in_perturbation_library=1"
            )
            _expressed_tfs = frozenset(r[0] for r in rows)
        except Exception:
            _expressed_tfs = frozenset()
    return _expressed_tfs


def _load_perturbed_tfs() -> frozenset:
    """TFs with at least one significant perturbation association in gsea_tf_table."""
    global _perturbed_tfs
    if _perturbed_tfs is None:
        try:
            rows = _query_db_raw(
                "SELECT DISTINCT gene_name FROM gsea_tf_table "
                "WHERE gene_set_collection='DE-hotspot_modules'"
            )
            _perturbed_tfs = frozenset(r[0] for r in rows)
        except Exception:
            _perturbed_tfs = frozenset()
    return _perturbed_tfs


def get_db():
    if "perturbseq_db" not in g:
        if not Path(DB_PATH).exists():
            abort(503)
        g.perturbseq_db = sqlite3.connect(DB_PATH)
        g.perturbseq_db.row_factory = sqlite3.Row
    return g.perturbseq_db


@perturbseq_bp.errorhandler(404)
def not_found(e):
    return render_template("perturbseq/404.html"), 404

@perturbseq_bp.errorhandler(500)
def internal_error(e):
    return render_template("perturbseq/500.html"), 500

@perturbseq_bp.errorhandler(503)
def db_unavailable(e):
    return render_template("perturbseq/db_error.html"), 503


@perturbseq_bp.context_processor
def inject_globals():
    return {'module_colors': MODULE_COLORS, 'last_updated': _fmt_time_ago(_LAST_COMMIT_TS)}


@perturbseq_bp.teardown_request
def close_db(exc):
    db = g.pop("perturbseq_db", None)
    if db is not None:
        db.close()


def _parse_sources_datasets(r) -> tuple:
    """Parse sources_str and datasets_str fields from a DB row."""
    sources = sorted(set((r['sources_str'] or '').split(','))) if r['sources_str'] else []
    seen_ds: set = set()
    datasets = []
    for entry in (r['datasets_str'] or '').split(','):
        parts = entry.split('|||')
        if len(parts) == 4 and parts[0] not in seen_ds:
            seen_ds.add(parts[0])
            datasets.append({'id': parts[0], 'name': parts[1], 'cell_type': parts[2], 'source': parts[3]})
    return sources, datasets


def rows_to_dicts(rows):
    out = []
    for r in rows:
        d = dict(r)
        for k, v in d.items():
            if isinstance(v, float):
                if math.isnan(v):
                    d[k] = None
                elif math.isinf(v):
                    d[k] = 1e308 if v > 0 else -1e308
        out.append(d)
    return out


_SPARK_TP_ORDER = _TP_ORDER

def _spark_pts(values, w=120, h=36, pad=3):
    """Return (line_d, area_d) smooth SVG path strings (Catmull-Rom spline) or (None, None)."""
    clean = [(i, v) for i, v in enumerate(values) if v is not None]
    if len(clean) < 2:
        return None, None
    vmin = min(v for _, v in clean)
    vmax = max(v for _, v in clean)
    vrange = vmax - vmin or 1
    n = len(values)
    bottom = h - pad
    pts = []
    for i, v in clean:
        x = pad + i / (n - 1) * (w - 2 * pad)
        y = bottom - (v - vmin) / vrange * (h - 2 * pad)
        pts.append((x, y))
    # Catmull-Rom → cubic bezier (clamped at endpoints)
    m = len(pts)
    d = f"M {pts[0][0]:.1f},{pts[0][1]:.1f}"
    for i in range(m - 1):
        p0 = pts[max(0, i - 1)]
        p1 = pts[i]
        p2 = pts[i + 1]
        p3 = pts[min(m - 1, i + 2)]
        cp1x = p1[0] + (p2[0] - p0[0]) / 6
        cp1y = p1[1] + (p2[1] - p0[1]) / 6
        cp2x = p2[0] - (p3[0] - p1[0]) / 6
        cp2y = p2[1] - (p3[1] - p1[1]) / 6
        d += f" C {cp1x:.1f},{cp1y:.1f} {cp2x:.1f},{cp2y:.1f} {p2[0]:.1f},{p2[1]:.1f}"
    area_d = d + f" L {pts[-1][0]:.1f},{bottom:.1f} L {pts[0][0]:.1f},{bottom:.1f} Z"
    return d, area_d


def _attach_sparklines(items, expr_map, key='name'):
    for item in items:
        tp_dict = expr_map.get(item[key], {})
        vals = [tp_dict.get(tp) for tp in _SPARK_TP_ORDER]
        item['spark_line'], item['spark_area'] = _spark_pts(vals)


def _all_modules_expr_map(db, source: str) -> dict:
    """Return {module_name: {timepoint: mean_tpm}} for all modules of given source.

    Replaces the pre-computed module_expression / submodule_expression tables.
    """
    rows = db.execute("""
        SELECT mt.module_name, be.sample_name AS sample, be.tpm
        FROM module_table mt
        JOIN gene_module_table gm ON mt.module_id = gm.module_id
        JOIN bulk_expression be ON gm.gene_id = be.gene_id
        WHERE mt.source = ? AND mt.module_name != 'unassigned' AND be.tpm IS NOT NULL
    """, (source,)).fetchall()
    by_mod_tp: dict = defaultdict(lambda: defaultdict(list))
    for mod, sample, tpm in rows:
        tp = _SAMPLE_TIMEPOINT_RE.sub('', sample)
        by_mod_tp[mod][tp].append(tpm)
    return {mod: {tp: sum(v) / len(v) for tp, v in by_tp.items()}
            for mod, by_tp in by_mod_tp.items()}


@perturbseq_bp.route("/")
def index():
    db = get_db()
    tfs, modules = _nav_lists(db)
    return render_template("perturbseq/index.html", tfs=tfs, modules=modules)


@perturbseq_bp.route("/all-modules")
def all_modules_page():
    tab = request.args.get("tab", "submodules")
    db = get_db()
    pub_counts = _load_gene_pub_counts()

    # ── Developmental Gene Clusters ──────────────────────────────────────────
    gc_rows = db.execute(
        "SELECT module_name, size FROM module_table "
        "WHERE source='mfuzz_k7' AND module_name != 'cluster_7' ORDER BY module_name"
    ).fetchall()
    gc_sizes = {_mfuzz_to_internal(r[0]): r[1] for r in gc_rows}

    gc_genes_raw = db.execute(
        "SELECT mt.module_name, gt.gene_name "
        "FROM gene_module_table gmt "
        "JOIN module_table mt ON gmt.module_id = mt.module_id "
        "JOIN gene_table gt ON gmt.gene_id = gt.gene_id "
        "WHERE mt.source = 'mfuzz_k7' AND mt.module_name != 'cluster_7'"
    ).fetchall()
    gc_genes_by_cluster: dict = defaultdict(list)
    for mod_name, gene in gc_genes_raw:
        gc_genes_by_cluster[_mfuzz_to_internal(mod_name)].append(gene)

    gc_regs_raw = db.execute(
        "SELECT module, gene_name, mean_NES, direction FROM gsea_tf_table "
        "WHERE module_collection='mfuzz_k7' "
        "AND gene_set_collection='ESC_DE-gene_clustering_data_var0.3_k7_top1000_DE' "
        "ORDER BY module, ABS(mean_NES) DESC"
    ).fetchall()
    gc_regs_by_cluster: dict = defaultdict(list)
    for mod_name, tf, nes, direction in gc_regs_raw:
        gc_regs_by_cluster[_mfuzz_to_internal(mod_name)].append({
            'tf': tf, 'nes': round(nes, 2),
            'direction': direction or ('up' if nes > 0 else 'down'),
        })

    clusters = []
    for cid in sorted(gc_sizes.keys()):
        genes = gc_genes_by_cluster[cid]
        notable = _top_genes_by_pub(genes, n=5)
        regs = gc_regs_by_cluster.get(cid, [])[:3]
        clusters.append({
            'id': cid,
            'display': _tc_display(cid),
            'size': gc_sizes[cid],
            'color': TC_CLUSTER_COLORS.get(cid, '#888'),
            'notable_genes': notable,
            'regulators': regs,
        })

    # ── Supermodules ──────────────────────────────────────────────────────────
    sup_rows = db.execute(
        "SELECT module_name, size FROM module_table "
        "WHERE source='hotspot_supermodule' AND module_name != 'unassigned' "
        "ORDER BY CAST(SUBSTR(module_name, 4) AS INTEGER)"
    ).fetchall()

    sup_genes_raw = db.execute(
        "SELECT gt.gene_name, mt.module_name "
        "FROM gene_module_table gm "
        "JOIN gene_table gt ON gm.gene_id = gt.gene_id "
        "JOIN module_table mt ON gm.module_id = mt.module_id "
        "WHERE mt.source='hotspot_supermodule' AND mt.module_name != 'unassigned'"
    ).fetchall()
    sup_genes_by_mod: dict = defaultdict(list)
    for gene, mod in sup_genes_raw:
        if gene:
            sup_genes_by_mod[mod].append(gene)

    sup_regs_raw = db.execute(
        "SELECT module, gene_name, mean_NES, direction FROM gsea_tf_table "
        "WHERE module_collection='hotspot_supermodule' "
        "AND gene_set_collection='DE-hotspot_modules' "
        "ORDER BY module, ABS(mean_NES) DESC"
    ).fetchall()
    sup_regs_by_mod: dict = defaultdict(list)
    for mod, tf, nes, direction in sup_regs_raw:
        sup_regs_by_mod[mod].append({
            'tf': tf, 'nes': round(nes, 2),
            'direction': direction or ('up' if nes > 0 else 'down'),
        })

    supermodules = []
    for mod_name, size in sup_rows:
        genes = sup_genes_by_mod.get(mod_name, [])
        notable = _top_genes_by_pub(genes, n=5)
        regs = sup_regs_by_mod.get(mod_name, [])[:3]
        supermodules.append({
            'name': mod_name,
            'size': size,
            'color': MODULE_COLORS.get(mod_name, '#888'),
            'notable_genes': notable,
            'regulators': regs,
        })

    # ── Submodules ────────────────────────────────────────────────────────────
    sub_rows = db.execute(
        "SELECT module_name, size, "
        "SUBSTR(module_name, 1, INSTR(module_name,'.')-1) AS supermodule "
        "FROM module_table WHERE source='hotspot_submodule' AND module_name != 'unassigned' "
        "ORDER BY CAST(SUBSTR(module_name, 4, INSTR(module_name,'.')-4) AS INTEGER), "
        "CAST(SUBSTR(module_name, INSTR(module_name,'.')+1) AS INTEGER)"
    ).fetchall()

    sub_genes_raw = db.execute(
        "SELECT gt.gene_name, mt.module_name "
        "FROM gene_module_table gm "
        "JOIN gene_table gt ON gm.gene_id = gt.gene_id "
        "JOIN module_table mt ON gm.module_id = mt.module_id "
        "WHERE mt.source='hotspot_submodule'"
    ).fetchall()
    sub_genes_by_mod: dict = defaultdict(list)
    for gene, mod in sub_genes_raw:
        if gene:
            sub_genes_by_mod[mod].append(gene)

    sub_regs_raw = db.execute(
        "SELECT module, gene_name, mean_NES, direction FROM gsea_tf_table "
        "WHERE module_collection='hotspot_submodule' "
        "AND gene_set_collection='DE-hotspot_modules' "
        "ORDER BY module, ABS(mean_NES) DESC"
    ).fetchall()
    sub_regs_by_mod: dict = defaultdict(list)
    for mod, tf, nes, direction in sub_regs_raw:
        sub_regs_by_mod[mod].append({
            'tf': tf, 'nes': round(nes, 2),
            'direction': direction or ('up' if nes > 0 else 'down'),
        })

    # ── GO enrichment terms (supermodules + submodules, bulk) ────────────────
    enr_raw = db.execute(
        "SELECT mt.module_name, e.term_name "
        "FROM go_module_enrichment e "
        "JOIN module_table mt ON e.module_id = mt.module_id "
        "WHERE mt.source IN ('hotspot_supermodule', 'hotspot_submodule') "
        "AND mt.module_name != 'unassigned' "
        "AND e.source IN ('GO:BP', 'GO:CC', 'GO:MF') "
        "ORDER BY mt.module_name, e.p_value"
    ).fetchall()
    enr_by_mod: dict = defaultdict(list)
    for mod, term in enr_raw:
        enr_by_mod[mod].append(term)

    submodules = []
    for mod_name, size, supermodule in sub_rows:
        genes = sub_genes_by_mod.get(mod_name, [])
        notable = _top_genes_by_pub(genes, n=5)
        regs = sub_regs_by_mod.get(mod_name, [])[:3]
        submodules.append({
            'name': mod_name,
            'supermodule': supermodule,
            'size': size,
            'color': MODULE_COLORS.get(supermodule, '#888'),
            'notable_genes': notable,
            'regulators': regs,
            'terms': enr_by_mod.get(mod_name, [])[:3],
        })

    for m in supermodules:
        m['terms'] = enr_by_mod.get(m['name'], [])[:3]

    # ── Expression sparklines (all three module types) ────────────────────────

    # Gene clusters — mean TPM per timepoint from bulk_expression
    gc_expr_map: dict = {
        _mfuzz_to_internal(k): v
        for k, v in _all_modules_expr_map(db, 'mfuzz_k7').items()
    }
    _attach_sparklines(clusters, gc_expr_map, key='id')

    # Supermodules — computed from bulk_expression
    sup_expr_map = _all_modules_expr_map(db, 'hotspot_supermodule')
    _attach_sparklines(supermodules, sup_expr_map)

    # Submodules — computed from bulk_expression
    sub_expr_map = _all_modules_expr_map(db, 'hotspot_submodule')
    _attach_sparklines(submodules, sub_expr_map)

    # ── Sort all three lists by developmental peak timing ────────────────────
    def _peak_tp_idx(tp_dict):
        pairs = [(i, tp_dict[tp]) for i, tp in enumerate(_SPARK_TP_ORDER)
                 if tp_dict.get(tp) is not None]
        return max(pairs, key=lambda x: x[1])[0] if pairs else len(_SPARK_TP_ORDER)

    for c in clusters:
        c['peak_idx'] = _peak_tp_idx(gc_expr_map.get(c['id'], {}))
    for m in supermodules:
        m['peak_idx'] = _peak_tp_idx(sup_expr_map.get(m['name'], {}))
    for m in submodules:
        m['peak_idx'] = _peak_tp_idx(sub_expr_map.get(m['name'], {}))

    clusters.sort(key=lambda c: c['peak_idx'])

    return render_template(
        "perturbseq/all_modules.html",
        tab=tab,
        clusters=clusters,
        supermodules=supermodules,
        submodules=submodules,
    )


@perturbseq_bp.route("/tf/<tf_name>")
def tf_page(tf_name):
    return redirect(url_for('perturbseq.gene_page', gene_name=tf_name))


@perturbseq_bp.route("/module/<name>")
def module_page(name):
    u = name.upper()
    if u[:2] == 'GC' and name[2:].isdigit():
        return _gc_page(name)
    if u[:2] == 'PC' and name[2:].isdigit():
        return _pc_page(name)
    if u[:3] == 'TC-' and name[3:].isdigit():
        return redirect(url_for('perturbseq.module_page', name='GC' + name[3:]), 301)
    if '.' in name:
        return _sub_page(name)
    return _super_page(name)


@perturbseq_bp.route("/submodule/<name>")
def submodule_page(name):
    return redirect(url_for('perturbseq.module_page', name=name), 301)


def _merge_pert_bind_edges(pert_rows, bind_map, id_key='tf'):
    """Merge perturbation + binding rows into a unified edge list with evidence field."""
    rows = []
    seen = set()
    for r in pert_rows:
        k = r[id_key]
        has_bind = k in bind_map
        seen.add(k)
        rows.append({
            id_key: k,
            'direction': 'Up' if r['mean_NES'] > 0 else 'Down',
            'mean_NES': r['mean_NES'], 'n_sig_gRNA': r['n_sig_gRNA'], 'padj': r['padj'],
            'binding': '✓' if has_bind else '',
            'odds_ratio': bind_map[k]['odds_ratio'] if has_bind else None,
            'binding_padj': bind_map[k].get('padj') if has_bind else None,
            'evidence': 'both' if has_bind else 'perturbation',
        })
    for k, b in bind_map.items():
        if k not in seen:
            rows.append({
                id_key: k, 'direction': '—', 'mean_NES': None,
                'n_sig_gRNA': None, 'padj': None, 'binding': '✓',
                'odds_ratio': b['odds_ratio'], 'binding_padj': b.get('padj'), 'evidence': 'binding',
            })
    return rows


def _build_gene_data_list(genes, pub_counts, lambert_tfs, perturbed_tfs, name_key='gene', pub_key='pub_count', empty_status=''):
    """Build gene data list with tf_status and publication count."""
    out = []
    for g in genes:
        if g in perturbed_tfs:
            status = 'Active TF'
        elif g in lambert_tfs:
            status = 'TF'
        else:
            status = empty_status
        out.append({name_key: g, 'tf_status': status, pub_key: pub_counts.get(g, 0)})
    return out


def _tf_dataset_ids(db, tf_gene_name: str) -> list:
    """Return dataset_ids for a TF — used to filter tf_peaks."""
    return [r[0] for r in db.execute(
        "SELECT dataset_id FROM tf_dataset_table WHERE tf_gene_name=?", (tf_gene_name,))]


def _bind_edges_for_tf(db, tf_gene_name: str, source_filter: str) -> dict:
    """Return {module_name: {odds_ratio, padj}} for a TF's binding enrichment.

    source_filter: 'hotspot_supermodule' or 'hotspot_submodule'
    """
    return {r['module']: r for r in rows_to_dicts(db.execute("""
        SELECT m.module_name AS module, MAX(e.odds_ratio) AS odds_ratio,
               MIN(e.padj_fisher) AS padj
        FROM tf_module_enrichment e
        JOIN module_table m ON e.module_id = m.module_id
        WHERE e.tf_gene_name = ? AND e.gene_set_collection = ?
        GROUP BY m.module_name
    """, (tf_gene_name, source_filter)))}


def _nav_lists(db):
    tfs = [r[0] for r in db.execute(
        "SELECT gene_name FROM gene_table WHERE in_perturbation_library=1 ORDER BY gene_name")]
    modules = [r[0] for r in db.execute(
        "SELECT module_name FROM module_table "
        "WHERE source='hotspot_supermodule' AND module_name != 'unassigned' ORDER BY module_name")]
    return tfs, modules


def _get_module_genes(db, module_id: int, ordered: bool = True) -> list:
    sql = ("SELECT gt.gene_name FROM gene_module_table gm "
           "JOIN gene_table gt ON gm.gene_id = gt.gene_id WHERE gm.module_id=?")
    if ordered:
        sql += " ORDER BY gt.gene_name"
    return [r[0] for r in db.execute(sql, (module_id,)).fetchall()]


def _get_module_description(db, module_id: int) -> dict | None:
    row = db.execute(
        "SELECT title, standard FROM module_description WHERE module_id=?",
        (module_id,)).fetchone()
    return dict(row) if row else None


def _get_module_go_highlights(db, module_id: int, limit: int = 5) -> list:
    return [r[0] for r in db.execute(
        "SELECT term_name FROM go_module_enrichment "
        "WHERE module_id=? AND source='GO:BP' ORDER BY p_value LIMIT ?",
        (module_id, limit)).fetchall()]


def _get_gene_membership(db, gene_id) -> dict:
    """Return {'submodule': str|None, 'supermodule': str|None}, normalising 'unassigned' to None."""
    rows = db.execute("""
        SELECT m.module_name, m.source
        FROM gene_module_table gm JOIN module_table m ON gm.module_id = m.module_id
        WHERE gm.gene_id=?
    """, (gene_id,)).fetchall()
    sub = next((r['module_name'] for r in rows if r['source'] == 'hotspot_submodule'), None)
    sup = next((r['module_name'] for r in rows if r['source'] == 'hotspot_supermodule'), None)
    if sub == 'unassigned': sub = None
    if sup == 'unassigned': sup = None
    if sub and not sup:
        sup = sub.rsplit('.', 1)[0] if '.' in sub else None
    return {'submodule': sub, 'supermodule': sup}


def _super_page(mod_name):
    db = get_db()
    mod_row = db.execute(
        "SELECT module_id FROM module_table WHERE module_name=? AND source='hotspot_supermodule'",
        (mod_name,)).fetchone()
    if not mod_row:
        return render_template("perturbseq/404.html", message=f"Module not found: {mod_name}"), 404
    module_id = mod_row['module_id']

    expr = _module_expr_by_timepoint(db, module_id)
    timepoints = [r['timepoint'] for r in expr]

    genes_list = _get_module_genes(db, module_id, ordered=False)
    expr_z = _compute_mean_zscore_profile(db, genes_list, timepoints)

    # Submodule list with top genes
    sub_rows = rows_to_dicts(db.execute(
        "SELECT module_id AS sub_module_id, module_name AS id, size AS n_genes "
        "FROM module_table WHERE source='hotspot_submodule' AND module_name LIKE ? ORDER BY size DESC",
        (mod_name + '.%',)))
    if sub_rows:
        pub_counts = _load_gene_pub_counts()
        sub_module_ids = [s['sub_module_id'] for s in sub_rows]
        ph = ','.join('?' * len(sub_module_ids))
        sub_gene_rows = db.execute(
            f"SELECT gm.module_id, gt.gene_name FROM gene_module_table gm "
            f"JOIN gene_table gt ON gm.gene_id = gt.gene_id "
            f"WHERE gm.module_id IN ({ph})",
            sub_module_ids).fetchall()
        sub_genes_map: dict = defaultdict(list)
        for mid, gname in sub_gene_rows:
            sub_genes_map[mid].append(gname)
        for s in sub_rows:
            genes = sub_genes_map.get(s['sub_module_id'], [])
            s['top_genes'] = _top_genes_by_pub(genes, n=5)
            s['color'] = MODULE_COLORS.get(mod_name, '#4a56b0')
            s['within_mean_z'] = 0.0
    submodules = sub_rows

    # TF perturbation edges
    pert_rows = rows_to_dicts(db.execute(
        "SELECT gene_name AS tf, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj "
        "FROM gsea_tf_table "
        "WHERE module=? AND module_collection='hotspot_supermodule' "
        "AND gene_set_collection='DE-hotspot_modules'",
        (mod_name,)))
    # TF binding edges: grouped by TF for this module
    bind_map = {}
    for r in db.execute("""
        SELECT e.tf_gene_name AS tf, MAX(e.odds_ratio) AS odds_ratio
        FROM tf_module_enrichment e
        WHERE e.module_id=? AND e.gene_set_collection='hotspot_supermodule'
        GROUP BY e.tf_gene_name
    """, (module_id,)).fetchall():
        bind_map[r['tf']] = {'tf': r['tf'], 'odds_ratio': r['odds_ratio']}
    top_tfs = _merge_pert_bind_edges(pert_rows, bind_map)
    for row in top_tfs:
        row['color'] = MODULE_COLORS.get(mod_name, '#4a56b0')
    top_tfs.sort(key=lambda r: abs(r['mean_NES'] or 0), reverse=True)

    notable = _top_genes_by_pub(genes_list)
    mod_desc_row = _get_module_description(db, module_id)
    highlight_terms = _get_module_go_highlights(db, module_id, limit=5)
    top_tf_names = [r['tf'] for r in top_tfs[:5] if r.get('tf')]
    desc = {
        'title': (mod_desc_row['title'] if mod_desc_row else None) or mod_name,
        'summary': (mod_desc_row['standard'] if mod_desc_row else None),
        'highlight_terms': highlight_terms,
        'highlight_genes': notable[:8],
        'top_tfs': top_tf_names,
        'similar_submodules': [],
    }
    # Tooltip data for submodule links in the table (description + expression sparkline)
    sub_module_ids = [s['sub_module_id'] for s in submodules] if submodules else []
    module_tooltips: dict = {}
    if sub_module_ids:
        ph = ','.join('?' * len(sub_module_ids))
        for _r in db.execute(
            f"SELECT mt.module_name, md.title, md.standard "
            f"FROM module_description md "
            f"JOIN module_table mt ON md.module_id = mt.module_id "
            f"WHERE mt.module_id IN ({ph})",
            sub_module_ids,
        ).fetchall():
            module_tooltips[_r['module_name']] = {
                'title': _r['title'],
                'desc':  _r['standard'],
            }
        sub_expr_map = _all_modules_expr_map(db, 'hotspot_submodule')
        for sub_name, tip in module_tooltips.items():
            tp_dict = sub_expr_map.get(sub_name, {})
            vals = [tp_dict.get(tp) for tp in _SPARK_TP_ORDER]
            line_d, area_d = _spark_pts(vals)
            if line_d:
                tip['spark_line'] = line_d
                tip['spark_area'] = area_d
    tfs, modules = _nav_lists(db)
    mod_color = MODULE_COLORS.get(mod_name, '#4a56b0')
    return render_template("perturbseq/module.html",
        module_type='super', name=mod_name, color=mod_color, card_color=mod_color,
        desc=desc, expr=expr, expr_z=expr_z, profile=None, parent=None, node=None,
        submodules=submodules, neighbors=[], genes=genes_list,
        notable_genes=notable,
        top_tfs=top_tfs, assoc=[], assoc_label='', enrichment=[],
        module_tooltips=module_tooltips,
        tfs=tfs, modules=modules)


def _sub_page(name):
    db = get_db()
    node_row = db.execute(
        "SELECT module_id, module_name AS id, size AS n_genes "
        "FROM module_table WHERE module_name=? AND source='hotspot_submodule'",
        (name,)).fetchone()
    if not node_row:
        return render_template("perturbseq/404.html", message=f"Submodule not found: {name}"), 404
    node = dict(node_row)
    module_id = node['module_id']
    supermodule = name.rsplit('.', 1)[0] if '.' in name else ''

    genes = _get_module_genes(db, module_id)
    pub_counts = _load_gene_pub_counts()
    lambert   = _load_lambert_tfs()
    perturbed = _load_perturbed_tfs()
    sub_gene_data = _build_gene_data_list(
        genes, pub_counts, lambert, perturbed,
        name_key='name', pub_key='pubs', empty_status='Gene')

    expr = _module_expr_by_timepoint(db, module_id)
    timepoints = [r['timepoint'] for r in expr]
    expr_z = _compute_mean_zscore_profile(db, genes, timepoints)

    highlight_terms = _get_module_go_highlights(db, module_id, limit=4)
    tf_rows = db.execute(
        "SELECT gene_name AS tf, direction FROM gsea_tf_table "
        "WHERE module=? AND module_collection='hotspot_submodule' "
        "AND gene_set_collection='DE-hotspot_modules' "
        "ORDER BY ABS(mean_NES) DESC LIMIT 5", (name,)).fetchall()
    top_tfs_auto = [r[0] for r in tf_rows]
    notable = _top_genes_by_pub(genes)

    summary_parts = []
    if highlight_terms:
        summary_parts.append(f"Enriched for {', '.join(highlight_terms[:3])}.")
    if tf_rows:
        tf_str = ', '.join(
            f"{tf} ({'↑' if direction == 'up' else '↓'})"
            for tf, direction in tf_rows[:3]
        )
        summary_parts.append(f"Top TF regulators: {tf_str}.")

    mod_desc_row = _get_module_description(db, module_id)
    desc = {
        'title': (mod_desc_row['title'] if mod_desc_row else None) or (highlight_terms[0] if highlight_terms else name),
        'summary': (mod_desc_row['standard'] if mod_desc_row else None) or (' '.join(summary_parts) or None),
        'highlight_terms': highlight_terms,
        'highlight_genes': notable,
        'top_tfs': top_tfs_auto,
        'similar_submodules': [],
    }
    tfs, modules = _nav_lists(db)
    card_color = MODULE_COLORS.get(supermodule, '#4a56b0')
    node['color'] = card_color
    node['supermodule'] = supermodule
    return render_template("perturbseq/module.html",
        module_type='sub', name=name, color=card_color, card_color=card_color,
        desc=desc, expr=expr, expr_z=expr_z, profile=None, parent=supermodule, node=node,
        submodules=[], neighbors=[], genes=genes,
        sub_gene_data=sub_gene_data,
        notable_genes=notable,
        top_tfs=top_tfs_auto, assoc=[], assoc_label='', enrichment=[],
        tfs=tfs, modules=modules)


def _gc_page(name):
    internal_id = 'gene_cluster_' + name[2:]
    if internal_id in HIDDEN_TC:
        return render_template("perturbseq/404.html", message=f"Cluster not found: {name}"), 404
    db = get_db()
    mfuzz_name = 'cluster_' + name[2:]
    tp_map = _all_modules_expr_map(db, 'mfuzz_k7').get(mfuzz_name, {})
    profile = [{'timepoint': tp, 'value': round(tp_map[tp], 5)}
               for tp in _SPARK_TP_ORDER if tp in tp_map]
    if not profile:
        return render_template("perturbseq/404.html", message=f"Cluster not found: {name}"), 404
    pert_rows = rows_to_dicts(db.execute("""
        SELECT gene_name AS tf, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj
        FROM gsea_tf_table
        WHERE module=? AND module_collection='mfuzz_k7'
        AND gene_set_collection='ESC_DE-gene_clustering_data_var0.3_k7_top1000_DE'
    """, (mfuzz_name,)).fetchall())
    mod_row = db.execute(
        "SELECT module_id FROM module_table WHERE module_name=? AND source='mfuzz_k7'",
        (mfuzz_name,)).fetchone()
    if mod_row:
        bind_map = {
            r[0]: {'odds_ratio': round(r[1], 3)}
            for r in db.execute("""
                SELECT e.tf_gene_name, MAX(e.odds_ratio) AS odds_ratio
                FROM tf_module_enrichment e
                WHERE e.module_id=? AND e.gene_set_collection='mfuzz' AND e.padj_fisher < 0.05
                GROUP BY e.tf_gene_name
            """, (mod_row['module_id'],)).fetchall()
        }
    else:
        bind_map = {}
    assoc = []
    top_tfs = _merge_pert_bind_edges(pert_rows, bind_map)
    top_tfs.sort(key=lambda r: abs(r['mean_NES'] or 0), reverse=True)
    genes = _get_module_genes(db, mod_row['module_id']) if mod_row else []
    pub_counts = _load_gene_pub_counts()
    lambert    = _load_lambert_tfs()
    perturbed  = _load_perturbed_tfs()
    gc_gene_data = _build_gene_data_list(
        genes, pub_counts, lambert, perturbed,
        name_key='name', pub_key='pubs', empty_status='Gene')
    timepoints = [r['timepoint'] for r in profile]
    expr_z = _compute_mean_zscore_profile(db, genes, timepoints)
    expr = _compute_mean_tpm_profile(db, genes, timepoints)
    desc = None
    if mod_row:
        gc_desc_row = _get_module_description(db, mod_row['module_id'])
        if gc_desc_row:
            highlight_terms = _get_module_go_highlights(db, mod_row['module_id'], limit=5)
            desc = {
                'title': gc_desc_row['title'] or name,
                'summary': gc_desc_row['standard'],
                'highlight_terms': highlight_terms,
                'highlight_genes': _top_genes_by_pub(genes)[:8],
                'top_tfs': [r['tf'] for r in top_tfs[:5] if r.get('tf')],
                'similar_submodules': [],
            }
    tfs, modules = _nav_lists(db)
    gc_color = TC_CLUSTER_COLORS.get(internal_id, '#4a56b0')
    return render_template("perturbseq/module.html",
        module_type='gc', name=name, color=gc_color, card_color=gc_color,
        desc=desc, expr=expr, expr_z=expr_z, profile=profile, parent=None, node=None,
        submodules=[], neighbors=[], genes=genes,
        gc_gene_data=gc_gene_data,
        notable_genes=_top_genes_by_pub(genes),
        top_tfs=top_tfs, assoc=assoc, assoc_label='Peak Clusters', enrichment=[],
        tfs=tfs, modules=modules)


def _pc_page(name):
    internal_id = 'peak_cluster_' + name[2:]
    if internal_id in HIDDEN_TC:
        return render_template("perturbseq/404.html", message=f"Cluster not found: {name}"), 404
    # Peak cluster data not available in consolidated DB
    return render_template("perturbseq/404.html",
                           message=f"Peak cluster pages are not available: {name}"), 404


@perturbseq_bp.route("/go/<path:term_name>")
def go_page(term_name):
    return _go_page(term_name)


def _go_page(term_name):
    db = get_db()

    # If given a bare GO accession (e.g. "GO:0007015"), resolve to go_name for redirect
    if re.match(r'^GO:\d+$', term_name):
        row = db.execute(
            "SELECT go_id, go_name FROM go_term_table WHERE go_accession=?", (term_name,)).fetchone()
        if row:
            return redirect(url_for('perturbseq.go_page',
                                    term_name=f"{row['go_name']} ({term_name})"), 301)
        return render_template("perturbseq/404.html",
                               message=f"GO term not found: {term_name}"), 404

    # Extract GO accession from term name string, e.g. "Cytoplasmic Translation (GO:0002181)"
    go_acc_match = re.search(r'\((GO:\d+)\)', term_name)
    go_accession = go_acc_match.group(1) if go_acc_match else None
    go_id = None  # integer PK in go_term_table

    if go_accession:
        go_row = db.execute(
            "SELECT go_id FROM go_term_table WHERE go_accession=? LIMIT 1", (go_accession,)).fetchone()
        if go_row:
            go_id = go_row['go_id']
    else:
        # Fallback: bare term name without accession
        go_row = db.execute(
            "SELECT go_id FROM go_term_table WHERE go_name=? LIMIT 1", (term_name,)).fetchone()
        if go_row:
            go_id = go_row['go_id']

    # Member genes via go_gene_table (go_id is integer PK)
    if go_id is not None:
        genes = [r[0] for r in db.execute("""
            SELECT g.gene_name FROM go_gene_table gg
            JOIN gene_table g ON gg.gene_id = g.gene_id
            WHERE gg.go_id=? ORDER BY g.gene_name
        """, (go_id,))]
    else:
        genes = []
    if not genes:
        return render_template("perturbseq/404.html",
                               message=f"GO term not found: {term_name}"), 404

    # Enriched modules via go_module_enrichment
    enriched_modules = []
    if go_id is not None:
        enriched_modules = rows_to_dicts(db.execute("""
            SELECT m.module_name AS module, m.source AS module_type,
                   e.p_value, e.term_size, e.intersection_size, e.query_size
            FROM go_module_enrichment e
            JOIN module_table m ON e.module_id = m.module_id
            WHERE e.go_id=? ORDER BY e.p_value
        """, (go_id,)))
        for row in enriched_modules:
            src = row['module_type']
            if src == 'hotspot_supermodule':
                row['module_type'] = 'supermodule'
                row['color'] = MODULE_COLORS.get(row['module'], '#4a56b0')
                row['url'] = url_for('perturbseq.module_page', name=row['module'])
            elif src == 'hotspot_submodule':
                row['module_type'] = 'submodule'
                sup = row['module'].rsplit('.', 1)[0] if '.' in row['module'] else ''
                row['color'] = MODULE_COLORS.get(sup, '#4a56b0')
                row['url'] = url_for('perturbseq.module_page', name=row['module'])
            else:
                row['color'] = '#4a56b0'
                row['url'] = '#'

    # TFs that perturb this GO term (via gsea_tf_table.go_id — integer FK)
    pert_rows = rows_to_dicts(db.execute(
        "SELECT gene_name AS tf, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj "
        "FROM gsea_tf_table WHERE go_id=? ORDER BY ABS(mean_NES) DESC",
        (go_id,))) if go_id is not None else []

    # TFs that bind this GO term gene set (via gmt_enrichment — gene_set uses accession string)
    bind_map = {r['tf']: r for r in rows_to_dicts(db.execute("""
        SELECT e.tf_gene_name AS tf, MAX(e.odds_ratio) AS odds_ratio
        FROM gmt_enrichment e
        WHERE e.gene_set LIKE ? AND e.gene_set_collection LIKE 'GO%'
        GROUP BY e.tf_gene_name
    """, (f'%({go_accession})%',)))} if go_accession else {}

    tf_rows = _merge_pert_bind_edges(pert_rows, bind_map)

    pub_counts = _load_gene_pub_counts()
    lambert = _load_lambert_tfs()
    perturbed = _load_perturbed_tfs()
    gene_data = _build_gene_data_list(genes, pub_counts, lambert, perturbed)
    gene_data.sort(key=lambda x: x['pub_count'], reverse=True)

    tfs, modules = _nav_lists(db)
    return render_template("perturbseq/go_term.html",
        term_name=term_name, go_id=go_id,
        genes=gene_data, n_genes=len(genes),
        enriched_modules=enriched_modules,
        tf_rows=tf_rows,
        tfs=tfs, modules=modules)


def _query_tf_elements(db, tf_name: str) -> list:
    """Return ATAC peaks bound by tf_name that also link to at least one gene."""
    ds_ids = _tf_dataset_ids(db, tf_name)
    if not ds_ids:
        return []
    ds_ph = ','.join('?' * len(ds_ids))
    rows = db.execute(f"""
        SELECT DISTINCT
            ap.atac_peak_id, ap.chr, ap.chrom_start AS start, ap.chrom_end AS end,
            gt.gene_name
        FROM tf_peaks          tp
        JOIN atac_tf_overlaps  ato ON tp.peak_id       = ato.peak_id
        JOIN atac_peak_table   ap  ON ato.atac_peak_id = ap.atac_peak_id
        JOIN (
            SELECT atac_peak_id, gene_id FROM atac_tss_links
            UNION ALL
            SELECT atac_peak_id, gene_id FROM multiome_atac_overlaps
        ) apg ON apg.atac_peak_id = ap.atac_peak_id
        JOIN gene_table        gt  ON apg.gene_id      = gt.gene_id
        WHERE tp.dataset_id IN ({ds_ph})
        ORDER BY ap.chr, ap.chrom_start, gt.gene_name
    """, ds_ids).fetchall()
    elements: dict = {}
    for r in rows:
        pk = r["atac_peak_id"]
        if pk not in elements:
            elements[pk] = {
                "atac_peak_id": pk,
                "chr":          _norm_chr(r["chr"]),
                "start":        r["start"],
                "end":          r["end"],
                "linked_genes": [],
            }
        if r["gene_name"] and r["gene_name"] not in elements[pk]["linked_genes"]:
            elements[pk]["linked_genes"].append(r["gene_name"])
    return list(elements.values())


@perturbseq_bp.route("/gene/<gene_name>")
def gene_page(gene_name):
    db = get_db()

    # Resolve gene_id via gene_table (or synonym fallback)
    gene_row = resolve_gene(db, gene_name)
    if not gene_row:
        return render_template("perturbseq/404.html",
                               message=f"Gene not found: {gene_name}"), 404
    gene_id = gene_row['gene_id']
    # Use canonical display name from DB in case query came via synonym
    gene_name = gene_row['gene_name']

    # Expression profile
    expr = _gene_expr_by_timepoint(db, gene_id)
    if not any(r['mean_tpm'] > 0 for r in expr):
        # gene exists in gene_table but has no expression data — still show the page
        pass

    # Module membership
    _mem = _get_gene_membership(db, gene_id)
    submodule = _mem['submodule']
    supermodule = _mem['supermodule']
    membership = _mem if (submodule or supermodule) else None

    # Developmental cluster membership from consolidated DB (mfuzz_k7 modules)
    tc_clusters = []
    for r in db.execute("""
        SELECT mt.module_name
        FROM gene_module_table gmt
        JOIN module_table mt ON gmt.module_id = mt.module_id
        JOIN gene_table gt ON gmt.gene_id = gt.gene_id
        WHERE mt.source = 'mfuzz_k7' AND gt.gene_name = ?
    """, (gene_name,)).fetchall():
        num = r[0].split('_')[1]          # "cluster_4" → "4"
        internal_id = 'gene_cluster_' + num
        tc_clusters.append({
            "cluster": internal_id,
            "label": 'GC' + num,
            "color": TC_CLUSTER_COLORS.get(internal_id, "#8090b0"),
            "url": url_for('perturbseq.module_page', name='GC' + num),
        })

    # TF status check
    is_tf = bool(db.execute(
        "SELECT 1 FROM gsea_tf_table WHERE gene_name=? AND gene_set_collection='DE-hotspot_modules' LIMIT 1",
        (gene_name,)).fetchone())

    # TF description from gene_table
    tf_desc_row = db.execute("SELECT description, Summary FROM gene_table WHERE gene_name=?", (gene_name,)).fetchone()
    gene_description = (tf_desc_row['description'] or None) if tf_desc_row else None
    tf_desc = {'name': gene_description, 'summary': tf_desc_row['Summary']} if (tf_desc_row and is_tf) else None
    gene_summary = (tf_desc_row['Summary'] if tf_desc_row and not is_tf else None)
    gene_aliases = [r['synonym'] for r in db.execute(
        "SELECT synonym FROM gene_synonym WHERE gene_id=? ORDER BY synonym_type, synonym",
        (gene_id,)).fetchall()]

    tf_reg_modules = []
    submodules_by_mod = {}
    pert_by_mod = {}
    tc_reg_clusters = []
    submodules_flat = []

    if is_tf:
        # Supermodule-level perturbation
        pert_rows = rows_to_dicts(db.execute(
            "SELECT module, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj "
            "FROM gsea_tf_table "
            "WHERE gene_name=? AND module_collection='hotspot_supermodule' "
            "AND gene_set_collection='DE-hotspot_modules'", (gene_name,)))
        # Supermodule-level binding
        bind_by_mod = _bind_edges_for_tf(db, gene_name, 'hotspot_supermodule')
        # Keep only modules with significant perturbation (|NES|>1, padj<0.05) or binding (OR>1, padj<0.05)
        pert_rows = [r for r in pert_rows if abs(r['mean_NES']) > 1 and r['padj'] < 0.05]
        bind_by_mod = {m: b for m, b in bind_by_mod.items() if b['odds_ratio'] > 1 and b['padj'] < 0.05}
        pert_mods = {r['module'] for r in pert_rows}
        pert_by_mod = {r['module']: r['mean_NES'] for r in pert_rows}

        tf_reg_modules = _merge_pert_bind_edges(pert_rows, bind_by_mod, id_key='module')
        for row in tf_reg_modules:
            row['color'] = MODULE_COLORS.get(row['module'], '#4a56b0')
        tf_reg_modules.sort(key=lambda r: abs(r['mean_NES'] or 0), reverse=True)

        # Submodule listing per regulated supermodule
        for mod in list(pert_mods | set(bind_by_mod.keys())):
            subs = rows_to_dicts(db.execute(
                "SELECT module_name AS id, size AS n_genes FROM module_table "
                "WHERE source='hotspot_submodule' AND module_name LIKE ? ORDER BY module_name",
                (mod + '.%',)))
            if subs:
                sup_color = MODULE_COLORS.get(mod, '#4a56b0')
                for s in subs:
                    s['supermodule_color'] = sup_color
                    s['color'] = sup_color
                    s['within_mean_z'] = 0.0
                submodules_by_mod[mod] = subs

        # Timecourse cluster regulation from consolidated DB
        HIDDEN = {'gene_cluster_7'}
        pert_gc = {
            _mfuzz_to_internal(r[0]): {'nes': r[1], 'padj': r[2], 'direction': r[3]}
            for r in db.execute("""
                SELECT module, mean_NES, min_padj, direction FROM gsea_tf_table
                WHERE gene_name=? AND module_collection='mfuzz_k7'
                AND gene_set_collection='ESC_DE-gene_clustering_data_var0.3_k7_top1000_DE'
            """, (gene_name,)).fetchall()
            if _mfuzz_to_internal(r[0]) not in HIDDEN
            and abs(r[1]) > 1 and r[2] < 0.05
        }
        # Binding enrichment from tf_module_enrichment (gene_set_collection='mfuzz')
        bind_gc = {
            _mfuzz_to_internal(mod): data
            for mod, data in _bind_edges_for_tf(db, gene_name, 'mfuzz').items()
            if _mfuzz_to_internal(mod) not in HIDDEN
            and data['odds_ratio'] > 1 and data['padj'] < 0.05
        }
        for gc in set(pert_gc.keys()) | set(bind_gc.keys()):
            p = pert_gc.get(gc)
            b = bind_gc.get(gc)
            evidence = 'both' if (p and b is not None) else ('perturbation' if p else 'binding')
            tc_reg_clusters.append({
                "cluster": gc,
                "label": gc.replace("gene_cluster_", "GC"),
                "nes": round(p['nes'], 3) if p else None,
                "padj": p['padj'] if p else None,
                "color": TC_CLUSTER_COLORS.get(gc, "#8090b0"),
                "url": url_for('perturbseq.module_page', name=_tc_display(gc)),
                "direction": ("Up" if p['nes'] > 0 else "Down") if p else "—",
                "binding": "✓" if b is not None else "",
                "odds_ratio": round(b['odds_ratio'], 3) if b else None,
                "binding_padj": b.get('padj') if b else None,
                "evidence": evidence,
            })
        tc_reg_clusters.sort(key=lambda r: abs(r['nes'] or 0), reverse=True)

        # Add expression sparkline paths for cluster node rendering
        if tc_reg_clusters:
            gc_expr = _all_modules_expr_map(db, 'mfuzz_k7')
            for entry in tc_reg_clusters:
                mfuzz_name = 'cluster_' + entry['cluster'].split('_')[2]
                tp_dict = gc_expr.get(mfuzz_name, {})
                vals = [tp_dict.get(tp) for tp in _SPARK_TP_ORDER]
                entry['spark_line'], entry['spark_area'] = _spark_pts(vals, w=88, h=28, pad=3)

        # Submodule-level perturbation + binding flat list
        sub_pert = {r['module']: r for r in rows_to_dicts(db.execute(
            "SELECT module, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj "
            "FROM gsea_tf_table WHERE gene_name=? AND module_collection='hotspot_submodule' "
            "AND gene_set_collection='DE-hotspot_modules'", (gene_name,)))
            if abs(r['mean_NES']) > 1 and r['padj'] < 0.05}
        sub_bind = {m: b for m, b in _bind_edges_for_tf(db, gene_name, 'hotspot_submodule').items()
                    if b['odds_ratio'] > 1 and b['padj'] < 0.05}
        all_sub_ids = set(sub_pert.keys()) | set(sub_bind.keys())
        if all_sub_ids:
            placeholders = ','.join('?' * len(all_sub_ids))
            sub_meta = {r['module_name']: r for r in rows_to_dicts(db.execute(
                f"SELECT module_name, size AS n_genes FROM module_table "
                f"WHERE source='hotspot_submodule' AND module_name IN ({placeholders})",
                list(all_sub_ids)))}
        else:
            sub_meta = {}
        for sub_id in sorted(all_sub_ids):
            meta = sub_meta.get(sub_id, {})
            sup = sub_id.rsplit('.', 1)[0] if '.' in sub_id else ''
            sup_color = MODULE_COLORS.get(sup, '#4a56b0')
            p = sub_pert.get(sub_id)
            b = sub_bind.get(sub_id)
            evidence = 'both' if (p and b) else ('perturbation' if p else 'binding')
            direction = ('Up' if p['mean_NES'] > 0 else 'Down') if p else None
            submodules_flat.append({
                'id': sub_id,
                'module': sup,
                'supermodule_color': sup_color,
                'n_genes': meta.get('n_genes', 0),
                'direction': direction,
                'mean_NES': p['mean_NES'] if p else None,
                'n_sig_gRNA': p['n_sig_gRNA'] if p else None,
                'padj': p['padj'] if p else None,
                'binding': '✓' if b else '',
                'odds_ratio': b['odds_ratio'] if b else None,
                'binding_padj': b.get('padj') if b else None,
                'evidence': evidence,
            })
        submodules_flat.sort(key=lambda x: abs(x['mean_NES'] or 0), reverse=True)

    card_color = MODULE_COLORS.get(supermodule or '', '#4a56b0')

    # Tooltip data for hotspot module pills
    _tip_mods = set()
    if submodule:  _tip_mods.add(submodule)
    if supermodule: _tip_mods.add(supermodule)
    for _m in tf_reg_modules:
        _tip_mods.add(_m['module'])
    module_tooltips: dict = {}
    if _tip_mods:
        _ph = ','.join('?' * len(_tip_mods))
        for _r in db.execute(
            f"SELECT mt.module_name, md.title, md.standard "
            f"FROM module_description md "
            f"JOIN module_table mt ON md.module_id = mt.module_id "
            f"WHERE mt.module_name IN ({_ph})",
            list(_tip_mods),
        ).fetchall():
            module_tooltips[_r['module_name']] = {
                'title': _r['title'],
                'desc': _r['standard'],
            }

    tfs, modules = _nav_lists(db)
    reg_elements     = _query_gene_elements(db, gene_name, gene_id=gene_id)
    tf_element_count = _count_tf_elements(db, gene_name) if is_tf else 0
    _locus_row = db.execute(
        'SELECT ap.chr, MIN(ap.chrom_start) AS locus_start, MAX(ap.chrom_end) AS locus_end '
        'FROM (SELECT atac_peak_id FROM atac_tss_links WHERE gene_id=? '
        '      UNION SELECT atac_peak_id FROM multiome_atac_overlaps WHERE gene_id=?) apg '
        'JOIN atac_peak_table ap ON apg.atac_peak_id = ap.atac_peak_id '
        'GROUP BY ap.chr ORDER BY (MAX(ap.chrom_end)-MIN(ap.chrom_start)) DESC LIMIT 1',
        (gene_id, gene_id)).fetchone()
    gene_locus = dict(_locus_row) if _locus_row else None
    if gene_locus:
        gene_locus['chr'] = _norm_chr(gene_locus['chr'])
        _tss_row = db.execute(
            "SELECT tss_position AS tss FROM gene_tss_table WHERE gene_id=? ORDER BY tss_id LIMIT 1",
            (gene_id,)).fetchone()
        gene_locus['tss'] = _tss_row['tss'] if _tss_row else None
    go_terms = [
        f"{r['go_name']} ({r['go_accession']})"
        for r in db.execute("""
            SELECT gt2.go_accession, gt2.go_name
            FROM go_gene_table gg
            JOIN gene_table g ON gg.gene_id = g.gene_id
            JOIN go_term_table gt2 ON gg.go_id = gt2.go_id
            WHERE g.gene_name = ? AND gt2.is_obsolete = 0
            ORDER BY gt2.go_name
        """, (gene_name,)).fetchall()
    ]

    return render_template("perturbseq/gene.html",
        gene=gene_name, expr=expr,
        membership=membership,
        tc_clusters=tc_clusters,
        is_tf=is_tf, tf_desc=tf_desc, gene_description=gene_description, gene_aliases=gene_aliases,
        card_color=card_color,
        tf_reg_modules=tf_reg_modules,
        submodules_by_mod=submodules_by_mod, pert_by_mod=pert_by_mod,
        submodules_flat=submodules_flat,
        tc_reg_clusters=tc_reg_clusters,
        tfs=tfs, modules=modules,
        reg_elements=reg_elements,
        tf_element_count=tf_element_count,
        gene_locus=gene_locus,
        gene_summary=gene_summary,
        go_terms=go_terms,
        module_tooltips=module_tooltips)


def _count_tf_elements(db, tf_name: str) -> int:
    # Pre-fetch dataset_ids to use idx_tf_peaks_dataset_id instead of scanning tf_dataset_table.
    # Returns exact count ≤501, or -1 as sentinel ("more than 501").
    ds_ids = _tf_dataset_ids(db, tf_name)
    if not ds_ids:
        return 0
    ph = ','.join('?' * len(ds_ids))
    rows = db.execute(f"""
        SELECT DISTINCT ato.atac_peak_id
        FROM tf_peaks         tp
        JOIN atac_tf_overlaps ato ON tp.peak_id = ato.peak_id
        WHERE tp.dataset_id IN ({ph})
          AND (EXISTS (SELECT 1 FROM atac_tss_links       WHERE atac_peak_id = ato.atac_peak_id)
            OR EXISTS (SELECT 1 FROM multiome_atac_overlaps WHERE atac_peak_id = ato.atac_peak_id))
        LIMIT 502
    """, ds_ids).fetchall()
    return len(rows) if len(rows) < 502 else -1


def _query_tf_elements_paged(db, tf_name: str, start: int, length: int,
                              q: str, order_col: int, order_dir: str) -> dict:
    """Paginated + filtered TF elements."""
    safe_dir = 'DESC' if order_dir.upper() == 'DESC' else 'ASC'
    order_map = {
        0: f'ap.chr {safe_dir}, ap.chrom_start {safe_dir}',
        1: f'(ap.chrom_end - ap.chrom_start) {safe_dir}',
    }
    order_by = order_map.get(order_col, f'ap.chr {safe_dir}, ap.chrom_start {safe_dir}')

    # Pre-fetch dataset_ids once — used for all subsequent queries
    ds_ids = _tf_dataset_ids(db, tf_name)
    if not ds_ids:
        return {'data': [], 'filtered': 0}
    ds_ph = ','.join('?' * len(ds_ids))

    if q:
        like = f'%{q}%'
        base_from = f"""
            FROM tf_peaks         tp
            JOIN atac_tf_overlaps ato ON tp.peak_id       = ato.peak_id
            JOIN atac_peak_table  ap  ON ato.atac_peak_id = ap.atac_peak_id
            JOIN (SELECT atac_peak_id, gene_id FROM atac_tss_links
                  UNION ALL
                  SELECT atac_peak_id, gene_id FROM multiome_atac_overlaps) apg
              ON apg.atac_peak_id = ap.atac_peak_id
            JOIN gene_table       gt  ON apg.gene_id      = gt.gene_id
            WHERE tp.dataset_id IN ({ds_ph})
            GROUP BY ap.atac_peak_id
            HAVING (ap.chr LIKE ? OR GROUP_CONCAT(DISTINCT gt.gene_name) LIKE ?)
        """
        params = ds_ids + [like, like]
        count_row = db.execute(
            f"SELECT COUNT(*) FROM (SELECT ap.atac_peak_id {base_from})", params
        ).fetchone()
        filtered_count = count_row[0] if count_row else 0
        rows = db.execute(f"""
            SELECT ap.atac_peak_id, ap.chr, ap.chrom_start AS start, ap.chrom_end AS end,
                   GROUP_CONCAT(DISTINCT gt.gene_name) AS genes_str
            {base_from}
            ORDER BY {order_by}
            LIMIT ? OFFSET ?
        """, params + [length, start]).fetchall()
        data = []
        for r in rows:
            genes = sorted(r['genes_str'].split(',')) if r['genes_str'] else []
            data.append({
                'atac_peak_id': r['atac_peak_id'],
                'chr':          _norm_chr(r['chr']),
                'start':        r['start'],
                'end':          r['end'],
                'linked_genes': genes,
            })
        return {'data': data, 'filtered': filtered_count}

    # Fast path (no search): use dataset_id IN (...) so idx_tf_peaks_dataset_id is hit.
    count_rows = db.execute(f"""
        SELECT DISTINCT ato.atac_peak_id
        FROM tf_peaks         tp
        JOIN atac_tf_overlaps ato ON tp.peak_id = ato.peak_id
        WHERE tp.dataset_id IN ({ds_ph})
          AND (EXISTS (SELECT 1 FROM atac_tss_links        WHERE atac_peak_id = ato.atac_peak_id)
            OR EXISTS (SELECT 1 FROM multiome_atac_overlaps WHERE atac_peak_id = ato.atac_peak_id))
        LIMIT 502
    """, ds_ids).fetchall()
    filtered_count = len(count_rows) if len(count_rows) < 502 else -1

    explicit_order = order_map.get(order_col) if order_col != 0 else None
    order_clause   = f"ORDER BY {explicit_order}" if explicit_order else ""

    page_rows = db.execute(f"""
        SELECT DISTINCT ap.atac_peak_id, ap.chr, ap.chrom_start AS start, ap.chrom_end AS end
        FROM tf_peaks         tp
        JOIN atac_tf_overlaps ato ON tp.peak_id       = ato.peak_id
        JOIN atac_peak_table  ap  ON ato.atac_peak_id = ap.atac_peak_id
        WHERE tp.dataset_id IN ({ds_ph})
          AND (EXISTS (SELECT 1 FROM atac_tss_links        WHERE atac_peak_id = ato.atac_peak_id)
            OR EXISTS (SELECT 1 FROM multiome_atac_overlaps WHERE atac_peak_id = ato.atac_peak_id))
        {order_clause}
        LIMIT ? OFFSET ?
    """, ds_ids + [length, start]).fetchall()

    if not page_rows:
        return {'data': [], 'filtered': filtered_count}

    peak_ids = [r['atac_peak_id'] for r in page_rows]
    placeholders = ','.join('?' * len(peak_ids))
    gene_rows = db.execute(f"""
        SELECT DISTINCT apg.atac_peak_id, gt.gene_name
        FROM (SELECT atac_peak_id, gene_id FROM atac_tss_links
              UNION ALL
              SELECT atac_peak_id, gene_id FROM multiome_atac_overlaps) apg
        JOIN gene_table gt ON apg.gene_id = gt.gene_id
        WHERE apg.atac_peak_id IN ({placeholders})
    """, peak_ids).fetchall()
    genes_by_peak: dict = {}
    for gr in gene_rows:
        genes_by_peak.setdefault(gr['atac_peak_id'], []).append(gr['gene_name'])

    data = []
    for r in page_rows:
        genes = sorted(genes_by_peak.get(r['atac_peak_id'], []))
        data.append({
            'atac_peak_id': r['atac_peak_id'],
            'chr':          _norm_chr(r['chr']),
            'start':        r['start'],
            'end':          r['end'],
            'linked_genes': genes,
        })
    return {'data': data, 'filtered': filtered_count}


@perturbseq_bp.route("/api/tf-elements/<tf_name>")
def api_tf_elements(tf_name):
    db = get_db()
    if 'start' in request.args:
        start      = int(request.args.get('start', 0))
        length     = int(request.args.get('length', 25))
        q          = request.args.get('q', '').strip()
        order_col  = int(request.args.get('order_col', 0))
        order_dir  = request.args.get('order_dir', 'asc')
        result = _query_tf_elements_paged(db, tf_name, start, length, q, order_col, order_dir)
        # When there's no search filter, total == filtered.  _count_tf_elements
        # is intentionally not called here; it can be very slow for genome-wide
        # binders (CTCF), and the filtered count already carries the -1 sentinel.
        total = result['filtered'] if not q else _count_tf_elements(db, tf_name)
        return jsonify({'data': result['data'], 'total': total, 'filtered': result['filtered']})
    return jsonify(_query_tf_elements(db, tf_name))


def _query_tf_linked_genes_paged(db, tf_name: str, start: int, length: int,
                                  q: str, order_col: int, order_dir: str,
                                  direction: str = 'all', source: str = '') -> dict:
    safe_dir = 'DESC' if order_dir.upper() == 'DESC' else 'ASC'
    order_map = {
        0: f'gt.gene_name {safe_dir}',
        1: f'sub.module_name {safe_dir}',
        2: f'n_datasets {safe_dir}',
        3: f'n_linked_peaks {safe_dir}',
        4: f'gdp.signed_pct_rank {safe_dir}',
        5: f'gdp.mean_coef {safe_dir}',
    }
    order_by = order_map.get(order_col, 'gdp.signed_pct_rank ASC')

    ds_ids = _tf_dataset_ids(db, tf_name)
    if not ds_ids:
        return {'data': [], 'filtered': 0, 'total': 0}
    ds_ph = ','.join('?' * len(ds_ids))

    # CTE that computes signed percentile rank (0=most downregulated, 1=most upregulated)
    # and abs percentile rank (for display bar). Only genes in the top or bottom 5% by
    # signed rank (i.e. strongest up- or down-regulated responders) are returned.
    de_cte = """
        WITH tf_grnas AS (
            SELECT grna_id FROM grna_table WHERE gene_name = ?
        ),
        gene_de AS (
            SELECT dr.gene_id, AVG(dr.coef) AS mean_coef
            FROM de_results dr
            WHERE dr.grna_id IN (SELECT grna_id FROM tf_grnas)
            GROUP BY dr.gene_id
        ),
        gene_de_pct AS (
            SELECT gene_id, mean_coef,
                   PERCENT_RANK() OVER (ORDER BY ABS(mean_coef)) AS de_pct_rank,
                   PERCENT_RANK() OVER (ORDER BY mean_coef)       AS signed_pct_rank
            FROM gene_de
        ),
        tf_bound_genes AS (
            SELECT DISTINCT apg.gene_id, tp.dataset_id
            FROM tf_peaks tp
            JOIN atac_tf_overlaps ato ON ato.peak_id = tp.peak_id
            JOIN (SELECT atac_peak_id, gene_id FROM atac_tss_links
                  UNION ALL
                  SELECT atac_peak_id, gene_id FROM multiome_atac_overlaps) apg
              ON apg.atac_peak_id = ato.atac_peak_id
    """

    if direction == 'up':
        de_filter = "gdp.signed_pct_rank >= 0.95"
    elif direction == 'down':
        de_filter = "gdp.signed_pct_rank <= 0.05"
    else:
        de_filter = "(gdp.signed_pct_rank >= 0.95 OR gdp.signed_pct_rank <= 0.05)"

    if source:
        source_clause = f"""AND tgl.gene_id IN (
            SELECT DISTINCT tbl.gene_id FROM tf_bound_genes tbl
            JOIN tf_dataset_table td2 ON tbl.dataset_id = td2.dataset_id
            WHERE td2.source = ?
        )"""
        src_params = [source]
    else:
        source_clause = ""
        src_params = []

    de_cte_closed = de_cte + f"            WHERE tp.dataset_id IN ({ds_ph})\n        )\n"

    sub_join = """
        LEFT JOIN (
            SELECT gm.gene_id, mo.module_name
            FROM gene_module_table gm
            JOIN module_table mo ON gm.module_id = mo.module_id
            WHERE mo.source = 'hotspot_submodule' AND mo.module_name != 'unassigned'
        ) sub ON tgl.gene_id = sub.gene_id
        LEFT JOIN (
            SELECT gm.gene_id, mo.module_name
            FROM gene_module_table gm
            JOIN module_table mo ON gm.module_id = mo.module_id
            WHERE mo.source = 'hotspot_supermodule' AND mo.module_name != 'unassigned'
        ) sup ON tgl.gene_id = sup.gene_id
    """

    total_row = db.execute(f"""
        {de_cte_closed}
        SELECT COUNT(DISTINCT tgl.gene_id)
        FROM tf_bound_genes tgl
        JOIN gene_de_pct gdp ON tgl.gene_id = gdp.gene_id
        WHERE {de_filter}
        {source_clause}
    """, [tf_name] + ds_ids + src_params).fetchone()
    total = total_row[0] if total_row else 0

    if q:
        like = f'%{q}%'
        having_clause = "HAVING (gt.gene_name LIKE ? OR sub.module_name LIKE ?)"
        q_params = [like, like]
        count_row = db.execute(f"""
            {de_cte_closed}
            SELECT COUNT(*) FROM (
                SELECT tgl.gene_id
                FROM tf_bound_genes tgl
                JOIN gene_table gt ON tgl.gene_id = gt.gene_id
                JOIN tf_dataset_table td ON tgl.dataset_id = td.dataset_id
                {sub_join}
                JOIN gene_de_pct gdp ON tgl.gene_id = gdp.gene_id
                WHERE {de_filter}
                {source_clause}
                GROUP BY tgl.gene_id, gt.gene_name
                {having_clause}
            )
        """, [tf_name] + ds_ids + src_params + q_params).fetchone()
        filtered = count_row[0] if count_row else 0
    else:
        having_clause = ""
        q_params = []
        filtered = total

    rows = db.execute(f"""
        {de_cte_closed}
        ,peak_counts AS (
            SELECT apg.gene_id, COUNT(DISTINCT ato.atac_peak_id) AS n_linked_peaks
            FROM tf_peaks tp
            JOIN atac_tf_overlaps ato ON ato.peak_id = tp.peak_id
            JOIN (SELECT atac_peak_id, gene_id FROM atac_tss_links
                  UNION ALL
                  SELECT atac_peak_id, gene_id FROM multiome_atac_overlaps) apg
              ON apg.atac_peak_id = ato.atac_peak_id
            WHERE tp.dataset_id IN ({ds_ph})
            GROUP BY apg.gene_id
        )
        SELECT gt.gene_name,
               sub.module_name AS submodule,
               sup.module_name AS supermodule,
               COUNT(DISTINCT tgl.dataset_id) AS n_datasets,
               GROUP_CONCAT(DISTINCT td.source) AS sources_str,
               GROUP_CONCAT(td.dataset_id || '|||' || td.dataset || '|||' || td.cell_type || '|||' || td.source) AS datasets_str,
               COALESCE(pc.n_linked_peaks, 0) AS n_linked_peaks,
               gdp.mean_coef,
               gdp.de_pct_rank,
               gdp.signed_pct_rank
        FROM tf_bound_genes tgl
        JOIN gene_table gt ON tgl.gene_id = gt.gene_id
        JOIN tf_dataset_table td ON tgl.dataset_id = td.dataset_id
        {sub_join}
        LEFT JOIN peak_counts pc ON tgl.gene_id = pc.gene_id
        JOIN gene_de_pct gdp ON tgl.gene_id = gdp.gene_id
        WHERE {de_filter}
        {source_clause}
        GROUP BY tgl.gene_id, gt.gene_name
        {having_clause}
        ORDER BY {order_by}
        LIMIT ? OFFSET ?
    """, [tf_name] + ds_ids + ds_ids + src_params + q_params + [length, start]).fetchall()

    data = []
    for r in rows:
        sources, datasets = _parse_sources_datasets(r)
        data.append({
            'gene_name':   r['gene_name'],
            'submodule':   r['submodule'],
            'supermodule': r['supermodule'],
            'n_datasets':  r['n_datasets'],
            'n_linked_peaks': r['n_linked_peaks'],
            'sources':     sources,
            'datasets':    datasets,
            'mean_coef':        r['mean_coef'],
            'de_pct_rank':      r['de_pct_rank'],
            'signed_pct_rank':  r['signed_pct_rank'],
        })
    return {'data': data, 'filtered': filtered, 'total': total}


@perturbseq_bp.route("/api/tf-linked-genes/<tf_name>")
def api_tf_linked_genes(tf_name):
    db = get_db()
    start     = int(request.args.get('start', 0))
    length    = int(request.args.get('length', 25))
    q         = request.args.get('q', '').strip()
    order_col = int(request.args.get('order_col', 4))
    order_dir = request.args.get('order_dir', 'asc')
    direction = request.args.get('direction', 'all')
    source    = request.args.get('source', '').strip()
    result = _query_tf_linked_genes_paged(db, tf_name, start, length, q, order_col, order_dir, direction, source)
    return jsonify(result)


@perturbseq_bp.route("/api/tf-linked-genes-all/<tf_name>")
def api_tf_linked_genes_all(tf_name):
    """Return all linked genes at once for client-side filtering/sorting."""
    db = get_db()
    ds_ids = _tf_dataset_ids(db, tf_name)
    if not ds_ids:
        return jsonify({'data': []})
    max_distance, link_types = _parse_link_filter()
    filter_sql, filter_params = _apg_filter_sql(max_distance, link_types)
    ds_ph = ','.join('?' * len(ds_ids))

    rows = db.execute(f"""
        WITH tf_grnas AS (
            SELECT grna_id FROM grna_table WHERE gene_name = ?
        ),
        gene_de AS (
            SELECT dr.gene_id, AVG(dr.coef) AS mean_coef, AVG(dr.z_coef) AS mean_z_coef
            FROM de_results dr
            WHERE dr.grna_id IN (SELECT grna_id FROM tf_grnas)
            GROUP BY dr.gene_id
        ),
        gene_de_pct AS (
            SELECT gene_id, mean_coef, mean_z_coef,
                   PERCENT_RANK() OVER (ORDER BY ABS(mean_coef)) AS de_pct_rank,
                   PERCENT_RANK() OVER (ORDER BY mean_coef)       AS signed_pct_rank
            FROM gene_de
        ),
        tf_bound_genes AS (
            SELECT DISTINCT apg.gene_id, tp.dataset_id
            FROM tf_peaks tp
            JOIN atac_tf_overlaps ato ON ato.peak_id = tp.peak_id
            JOIN (SELECT atac_peak_id, gene_id, distance_bp, NULL AS mao_lt FROM atac_tss_links
                  UNION ALL
                  SELECT atac_peak_id, gene_id, distance_to_tss, link_type FROM multiome_atac_overlaps) apg
              ON apg.atac_peak_id = ato.atac_peak_id
            WHERE tp.dataset_id IN ({ds_ph}){filter_sql}
        ),
        peak_type_counts AS (
            SELECT apg.gene_id,
                   COUNT(DISTINCT ato.atac_peak_id) AS n_peaks,
                   COUNT(DISTINCT CASE WHEN apg.distance_bp <= 1000 AND apg.mao_lt IS NULL
                                       THEN ato.atac_peak_id END) AS n_proximal,
                   COUNT(DISTINCT CASE WHEN apg.distance_bp > 1000 AND apg.mao_lt IS NULL
                                       THEN ato.atac_peak_id END) AS n_distal,
                   COUNT(DISTINCT CASE WHEN apg.mao_lt IS NOT NULL
                                       THEN ato.atac_peak_id END) AS n_multiome
            FROM tf_peaks tp
            JOIN atac_tf_overlaps ato ON ato.peak_id = tp.peak_id
            JOIN (SELECT atac_peak_id, gene_id, distance_bp, NULL AS mao_lt FROM atac_tss_links
                  UNION ALL
                  SELECT atac_peak_id, gene_id, distance_to_tss, link_type FROM multiome_atac_overlaps) apg
              ON apg.atac_peak_id = ato.atac_peak_id
            WHERE tp.dataset_id IN ({ds_ph}){filter_sql}
            GROUP BY apg.gene_id
        )
        SELECT gt.gene_name,
               sub.module_name AS submodule,
               sup.module_name AS supermodule,
               COUNT(DISTINCT tgl.dataset_id) AS n_datasets,
               GROUP_CONCAT(DISTINCT td.source) AS sources_str,
               GROUP_CONCAT(td.dataset_id || '|||' || td.dataset || '|||' || td.cell_type || '|||' || td.source) AS datasets_str,
               COALESCE(ptc.n_peaks,        0)   AS n_peaks,
               COALESCE(ptc.n_proximal,     0)   AS n_proximal,
               COALESCE(ptc.n_distal,       0)   AS n_distal,
               COALESCE(ptc.n_multiome,     0)   AS n_multiome,
               gdp.mean_coef,
               gdp.mean_z_coef,
               gdp.de_pct_rank,
               gdp.signed_pct_rank
        FROM tf_bound_genes tgl
        JOIN gene_table gt ON tgl.gene_id = gt.gene_id
        JOIN tf_dataset_table td ON tgl.dataset_id = td.dataset_id
        LEFT JOIN (
            SELECT gm.gene_id, mo.module_name FROM gene_module_table gm
            JOIN module_table mo ON gm.module_id = mo.module_id
            WHERE mo.source = 'hotspot_submodule' AND mo.module_name != 'unassigned'
        ) sub ON tgl.gene_id = sub.gene_id
        LEFT JOIN (
            SELECT gm.gene_id, mo.module_name FROM gene_module_table gm
            JOIN module_table mo ON gm.module_id = mo.module_id
            WHERE mo.source = 'hotspot_supermodule' AND mo.module_name != 'unassigned'
        ) sup ON tgl.gene_id = sup.gene_id
        LEFT JOIN peak_type_counts ptc ON tgl.gene_id = ptc.gene_id
        JOIN gene_de_pct gdp ON tgl.gene_id = gdp.gene_id
        WHERE (gdp.signed_pct_rank >= 0.95 OR gdp.signed_pct_rank <= 0.05)
        GROUP BY tgl.gene_id, gt.gene_name
        ORDER BY gdp.signed_pct_rank ASC
    """, [tf_name] + ds_ids + filter_params + ds_ids + filter_params).fetchall()

    data = []
    for r in rows:
        sources, datasets = _parse_sources_datasets(r)
        data.append({
            'gene_name':        r['gene_name'],
            'submodule':        r['submodule'],
            'supermodule':      r['supermodule'],
            'n_datasets':       r['n_datasets'],
            'n_peaks':          r['n_peaks'],
            'n_proximal':       r['n_proximal'],
            'n_distal':         r['n_distal'],
            'n_multiome':       r['n_multiome'],
            'sources':          sources,
            'datasets':         datasets,
            'mean_coef':        r['mean_coef'],
            'mean_z_coef':      r['mean_z_coef'],
            'de_pct_rank':      r['de_pct_rank'],
            'signed_pct_rank':  r['signed_pct_rank'],
        })
    return jsonify({'data': data})


_gene_coef_cache: dict = {}   # tf_name → {gene_id: mean_coef}
_coef_rug_cache:  dict = {}   # tf_name → linked_by_type payload


def _get_gene_coef(db, tf_name: str) -> dict:
    """Return {gene_id: mean_coef (2dp)} for all genes measured in the TF KD experiment."""
    if tf_name in _gene_coef_cache:
        return _gene_coef_cache[tf_name]
    coef_rows = db.execute("""
        SELECT dr.gene_id, AVG(dr.coef) AS mean_coef
        FROM de_results dr
        JOIN grna_table gr ON dr.grna_id = gr.grna_id
        WHERE gr.gene_name = ?
        GROUP BY dr.gene_id
    """, (tf_name,)).fetchall()
    gene_coef: dict = {r['gene_id']: round(r['mean_coef'], 2)
                       for r in coef_rows if r['mean_coef'] is not None}
    _gene_coef_cache[tf_name] = gene_coef
    return gene_coef


def _get_linked_by_type(db, tf_name: str, gene_coef: dict) -> dict:
    """Return linked_by_type payload for rug plot (only top/bottom 5% DE responders)."""
    if tf_name in _coef_rug_cache:
        return _coef_rug_cache[tf_name]

    empty = {'proximity': [], 'multiome': [], 'multiome_hicar': []}
    if not gene_coef:
        return empty

    coefs_sorted = sorted(gene_coef.items(), key=lambda kv: kv[1])
    n = len(coefs_sorted)
    lo_ids = {gid for gid, _ in coefs_sorted[:max(1, int(n * 0.05))]}
    hi_ids = {gid for gid, _ in coefs_sorted[max(0, int(n * 0.95)):]}
    extreme_ids = list(lo_ids | hi_ids)
    if not extreme_ids:
        _coef_rug_cache[tf_name] = empty
        return empty
    ex_ph = ','.join('?' * len(extreme_ids))

    linked_rows = db.execute(f"""
        SELECT DISTINCT
            apg.gene_id,
            gt.gene_name,
            CASE
                WHEN apg.mao_lt IN ('HiCAR+TSS', 'HiCAR') THEN 'multiome_hicar'
                WHEN apg.mao_lt = 'TSS'                    THEN 'multiome'
                ELSE 'proximity'
            END AS evidence_type
        FROM tf_dataset_table td
        JOIN tf_peaks         tp  ON tp.dataset_id   = td.dataset_id
        JOIN atac_tf_overlaps ato ON ato.peak_id     = tp.peak_id
        JOIN (SELECT atac_peak_id, gene_id, distance_bp, NULL AS mao_lt FROM atac_tss_links
              UNION ALL
              SELECT atac_peak_id, gene_id, distance_to_tss, link_type FROM multiome_atac_overlaps) apg
          ON apg.atac_peak_id = ato.atac_peak_id
        JOIN gene_table       gt  ON gt.gene_id      = apg.gene_id
        WHERE td.tf_gene_name = ?
        AND apg.gene_id IN ({ex_ph})
    """, (tf_name, *extreme_ids)).fetchall()

    buckets: dict = {'proximity': {}, 'multiome': {}, 'multiome_hicar': {}}
    for r in linked_rows:
        coef = gene_coef.get(r['gene_id'])
        if coef is not None:
            lt = r['evidence_type']
            bucket = round(coef, 2)
            if bucket not in buckets[lt]:
                buckets[lt][bucket] = []
            buckets[lt][bucket].append(r['gene_name'])

    _MAX_GENES = 10
    linked_by_type = {
        k: [{'coef': c, 'genes': sorted(gs)[:_MAX_GENES], 'n': len(gs)}
            for c, gs in sorted(v.items())]
        for k, v in buckets.items()
    }
    _coef_rug_cache[tf_name] = linked_by_type
    return linked_by_type


@perturbseq_bp.route("/api/tf/<tf_name>/coef-density")
def api_tf_coef_density(tf_name):
    """Fast endpoint — returns all-gene mean_coefs for KDE (~2s even for large TFs)."""
    db = get_db()
    ds_ids = _tf_dataset_ids(db, tf_name)
    if not ds_ids:
        return jsonify({'all_coefs': [], 'bound_coefs': []})
    gene_coef = _get_gene_coef(db, tf_name)
    max_distance, link_types = _parse_link_filter()
    filter_sql, filter_params = _apg_filter_sql(max_distance, link_types)
    ds_ph = ','.join('?' * len(ds_ids))
    bound_rows = db.execute(f"""
        SELECT DISTINCT apg.gene_id, gt.gene_name
        FROM tf_peaks tp
        JOIN atac_tf_overlaps ato ON ato.peak_id = tp.peak_id
        JOIN (SELECT atac_peak_id, gene_id, distance_bp, NULL AS mao_lt FROM atac_tss_links
              UNION ALL
              SELECT atac_peak_id, gene_id, distance_to_tss, link_type FROM multiome_atac_overlaps) apg
          ON apg.atac_peak_id = ato.atac_peak_id
        JOIN gene_table gt ON apg.gene_id = gt.gene_id
        WHERE tp.dataset_id IN ({ds_ph}){filter_sql}
    """, ds_ids + filter_params).fetchall()
    bound_id_to_name = {r['gene_id']: r['gene_name'] for r in bound_rows}
    bound_genes = [{'gene_name': bound_id_to_name[gid], 'coef': coef}
                   for gid, coef in gene_coef.items()
                   if gid in bound_id_to_name]
    return jsonify({'all_coefs': list(gene_coef.values()), 'bound_genes': bound_genes})


@perturbseq_bp.route("/api/tf/<tf_name>/coef-density/rug")
def api_tf_coef_density_rug(tf_name):
    """Slow endpoint — returns rug data (linked genes by type).  May take ~50s first call."""
    db = get_db()
    if not _tf_dataset_ids(db, tf_name):
        return jsonify({'linked_by_type': {}})
    gene_coef = _get_gene_coef(db, tf_name)
    return jsonify({'linked_by_type': _get_linked_by_type(db, tf_name, gene_coef)})


@perturbseq_bp.route("/api/gene/<gene_name>")
def api_gene(gene_name):
    db = get_db()
    gene_row = resolve_gene(db, gene_name)
    if not gene_row:
        return jsonify({"error": "not found"}), 404
    expr = _gene_expr_by_timepoint(db, gene_row['gene_id'])
    membership = _get_gene_membership(db, gene_row['gene_id'])
    return jsonify({"gene": gene_name, "expr": expr, "membership": membership})


@perturbseq_bp.route("/api/submodule/<name>")
def api_submodule(name):
    db = get_db()
    node_row = db.execute(
        "SELECT module_id, module_name AS id, size AS n_genes FROM module_table "
        "WHERE module_name=? AND source='hotspot_submodule'", (name,)).fetchone()
    if not node_row:
        return jsonify({"error": "not found"}), 404
    node = dict(node_row)
    node['color'] = MODULE_COLORS.get(name.rsplit('.', 1)[0] if '.' in name else '', '#4a56b0')
    node['within_mean_z'] = 0.0
    genes = _get_module_genes(db, node['module_id'])
    return jsonify({"node": node, "genes": genes, "neighbors": []})


@perturbseq_bp.route("/api/edges")
def api_edges():
    db = get_db()
    min_strength = float(request.args.get("min_strength", 0))
    evidence = request.args.get("evidence", "all")
    tf_filter = request.args.get("tf", "").upper()
    mod_filter = request.args.get("module", "").upper()

    expressed = _load_expressed_tfs()

    pert = rows_to_dicts(db.execute(
        "SELECT gene_name AS tf, module, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj "
        "FROM gsea_tf_table WHERE gene_set_collection='DE-hotspot_modules'").fetchall())
    bind_map = {}
    for r in db.execute("""
        SELECT e.tf_gene_name AS tf, m.module_name AS module,
               MAX(e.odds_ratio) AS odds_ratio, MIN(e.padj_fisher) AS padj
        FROM tf_module_enrichment e
        JOIN module_table m ON e.module_id = m.module_id
        WHERE e.gene_set_collection IN ('hotspot_supermodule', 'hotspot_submodule')
        GROUP BY e.tf_gene_name, m.module_name
    """).fetchall():
        if expressed and r["tf"] not in expressed:
            continue
        bind_map[(r["tf"], r["module"])] = {"odds_ratio": r["odds_ratio"], "padj": r["padj"]}

    edges = []
    seen = set()
    for e in pert:
        key = (e["tf"], e["module"])
        seen.add(key)
        has_bind = key in bind_map
        strength = abs(e["mean_NES"] or 0)
        if has_bind:
            strength = max(strength, bind_map[key]["odds_ratio"] or 0)
        if strength < min_strength:
            continue
        if evidence == "both" and not has_bind:
            continue
        if evidence == "binding" and not has_bind:
            continue
        if tf_filter and tf_filter not in e["tf"].upper():
            continue
        if mod_filter and mod_filter not in e["module"].upper():
            continue
        edges.append({
            "tf": e["tf"], "module": e["module"],
            "perturbation": True, "binding": has_bind,
            "mean_NES": e["mean_NES"], "strength": strength,
        })

    if evidence != "perturbation":
        for (tf, mod), b in bind_map.items():
            if (tf, mod) in seen:
                continue
            strength = b["odds_ratio"] or 0
            if strength < min_strength:
                continue
            if evidence == "both":
                continue
            if tf_filter and tf_filter not in tf.upper():
                continue
            if mod_filter and mod_filter not in mod.upper():
                continue
            edges.append({
                "tf": tf, "module": mod,
                "perturbation": False, "binding": True,
                "mean_NES": None, "strength": strength,
            })

    return jsonify(edges)


@perturbseq_bp.route("/api/tf/<tf_name>")
def api_tf(tf_name):
    db = get_db()
    pert = rows_to_dicts(db.execute(
        "SELECT module, module_collection, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj, direction "
        "FROM gsea_tf_table WHERE gene_name=? AND gene_set_collection='DE-hotspot_modules'",
        (tf_name,)).fetchall())
    bind = rows_to_dicts(db.execute("""
        SELECT m.module_name AS module, m.source AS module_collection,
               MAX(e.odds_ratio) AS odds_ratio, MIN(e.padj_fisher) AS padj
        FROM tf_module_enrichment e
        JOIN module_table m ON e.module_id = m.module_id
        WHERE e.tf_gene_name = ? AND e.gene_set_collection IN ('hotspot_supermodule','hotspot_submodule')
        GROUP BY m.module_name
    """, (tf_name,)).fetchall())
    return jsonify({"perturbation": pert, "binding": bind})


@perturbseq_bp.route("/api/module/<mod_name>")
def api_module(mod_name):
    db = get_db()
    mod_row = db.execute(
        "SELECT module_id FROM module_table WHERE module_name=?", (mod_name,)).fetchone()
    if not mod_row:
        return jsonify({"perturbation": [], "binding": []})
    module_id = mod_row['module_id']
    pert = rows_to_dicts(db.execute(
        "SELECT gene_name AS tf, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj, direction "
        "FROM gsea_tf_table WHERE module=? AND gene_set_collection='DE-hotspot_modules'",
        (mod_name,)).fetchall())
    bind = rows_to_dicts(db.execute("""
        SELECT e.tf_gene_name AS tf, MAX(e.odds_ratio) AS odds_ratio, MIN(e.padj_fisher) AS padj
        FROM tf_module_enrichment e
        WHERE e.module_id=? AND e.gene_set_collection IN ('hotspot_supermodule','hotspot_submodule')
        GROUP BY e.tf_gene_name
    """, (module_id,)).fetchall())
    return jsonify({"perturbation": pert, "binding": bind})


@perturbseq_bp.route("/api/module/<mod_name>/detail")
def api_module_detail(mod_name):
    db = get_db()
    mod_row = db.execute(
        "SELECT module_id FROM module_table WHERE module_name=?", (mod_name,)).fetchone()
    if not mod_row:
        return jsonify({"enrichment": [], "genes": [], "pub_counts": {},
                        "gene_tf_status": {}, "gene_submodule": {}})
    module_id = mod_row['module_id']
    enrichment = rows_to_dicts(db.execute(
        "SELECT term_id, source, term_name, p_value, term_size, intersection_size "
        "FROM go_module_enrichment WHERE module_id=? AND source IN ('GO:BP','GO:CC','GO:MF')",
        (module_id,)).fetchall())
    genes = _get_module_genes(db, module_id)
    pub_counts = _load_gene_pub_counts()
    gene_pubs = {g: pub_counts.get(g, 0) for g in genes}
    lambert   = _load_lambert_tfs()
    perturbed = _load_perturbed_tfs()
    gene_tf_status = {
        row['gene']: row['tf_status']
        for row in _build_gene_data_list(genes, pub_counts, lambert, perturbed, empty_status='Gene')
    }
    # Map genes to their submodule (if any)
    sub_rows = db.execute(
        "SELECT module_id AS sub_id, module_name FROM module_table "
        "WHERE source='hotspot_submodule' AND module_name LIKE ?",
        (mod_name + '.%',)).fetchall()
    sub_color = MODULE_COLORS.get(mod_name, '#4a56b0')
    gene_submodule = {}
    for sub_row in sub_rows:
        for gr in db.execute(
            "SELECT gt.gene_name FROM gene_module_table gm "
            "JOIN gene_table gt ON gm.gene_id = gt.gene_id WHERE gm.module_id=?",
                (sub_row['sub_id'],)).fetchall():
            gene_submodule[gr[0]] = {"id": sub_row['module_name'], "color": sub_color}
    return jsonify({"enrichment": enrichment, "genes": genes,
                    "pub_counts": gene_pubs, "gene_tf_status": gene_tf_status,
                    "gene_submodule": gene_submodule})


@perturbseq_bp.route("/api/module/<mod_name>/genes")
def api_module_genes(mod_name):
    db = get_db()
    mod_row = db.execute(
        "SELECT module_id FROM module_table WHERE module_name=?", (mod_name,)).fetchone()
    if not mod_row:
        return jsonify([])
    genes = _get_module_genes(db, mod_row['module_id'])
    return jsonify(genes)


@perturbseq_bp.route("/api/module/<mod_name>/tooltip")
def api_module_tooltip(mod_name):
    db = get_db()
    # GC display names (GC3) map to mfuzz_k7 module names (cluster_3)
    lookup = ('cluster_' + mod_name[2:]) if re.match(r'^GC\d+$', mod_name) else mod_name
    mod_row = db.execute(
        "SELECT module_id FROM module_table WHERE module_name=?", (lookup,)).fetchone()
    if not mod_row:
        return jsonify({'title': None, 'desc': None})
    desc_row = _get_module_description(db, mod_row['module_id'])
    result = {
        'title': desc_row['title'] if desc_row else None,
        'desc':  desc_row['standard'] if desc_row else None,
    }
    if re.match(r'^GC\d+$', mod_name):
        internal_id = 'gene_cluster_' + mod_name[2:]
        genes = [r[0] for r in db.execute(
            "SELECT gt.gene_name FROM gene_module_table gmt "
            "JOIN module_table mt ON gmt.module_id = mt.module_id "
            "JOIN gene_table gt ON gmt.gene_id = gt.gene_id "
            "WHERE mt.module_name=? AND mt.source='mfuzz_k7'",
            (lookup,)).fetchall()]
        result['expr']  = _compute_mean_zscore_profile(db, genes, _TP_ORDER)
        result['color'] = TC_CLUSTER_COLORS.get(internal_id, '#4a56b0')
    return jsonify(result)


@perturbseq_bp.route("/api/gene/<gene_name>/tooltip")
def api_gene_tooltip(gene_name):
    db = get_db()
    gene_row = resolve_gene(db, gene_name)
    if not gene_row:
        return jsonify({'full_name': None, 'summary': None, 'submodule': None,
                        'supermodule': None, 'is_tf': False})

    desc_row = db.execute(
        "SELECT description, Summary FROM gene_table WHERE gene_id=?",
        (gene_row['gene_id'],)
    ).fetchone()
    gene_name_full = (desc_row['description'] or None) if desc_row else None
    summary        = (desc_row['Summary']     or None) if desc_row else None

    _mem = _get_gene_membership(db, gene_row['gene_id'])
    submodule   = _mem['submodule']
    supermodule = _mem['supermodule']

    is_tf = bool(db.execute(
        "SELECT 1 FROM gsea_tf_table WHERE gene_name=? AND gene_set_collection='DE-hotspot_modules' LIMIT 1",
        (gene_name,)
    ).fetchone())

    return jsonify({
        'full_name':   gene_name_full,
        'summary':     summary,
        'submodule':   submodule,
        'supermodule': supermodule,
        'is_tf':       is_tf,
    })


@perturbseq_bp.route("/api/gene/<gene_name>/grna-coefs")
def api_gene_grna_coefs(gene_name):
    db = get_db()
    gene_row = resolve_gene(db, gene_name)
    if not gene_row:
        return jsonify({'coefs': []})

    # Precompute TFs that have binding evidence for this gene by starting from
    # the gene's linked ATAC peaks (small set) rather than scanning per-row.
    rows = db.execute("""
        WITH gene_peaks AS (
            SELECT atac_peak_id FROM atac_tss_links WHERE gene_id = ?
            UNION
            SELECT atac_peak_id FROM multiome_atac_overlaps WHERE gene_id = ?
        ),
        binding AS (
            SELECT DISTINCT tdt.tf_gene_name
            FROM gene_peaks gp
            JOIN atac_tf_overlaps ato ON ato.atac_peak_id = gp.atac_peak_id
            JOIN tf_peaks tp ON tp.peak_id = ato.peak_id
            JOIN tf_dataset_table tdt ON tdt.dataset_id = tp.dataset_id
        )
        SELECT gr.gene_name AS tf_name, gr.grna_name, gr.active, dr.coef,
               PERCENT_RANK() OVER (ORDER BY dr.coef) AS signed_pct_rank,
               CASE WHEN b.tf_gene_name IS NOT NULL THEN 1 ELSE 0 END AS has_binding
        FROM de_results dr
        JOIN grna_table gr ON dr.grna_id = gr.grna_id
        LEFT JOIN binding b ON b.tf_gene_name = gr.gene_name
        WHERE dr.gene_id = ?
        ORDER BY dr.coef
    """, (gene_row['gene_id'], gene_row['gene_id'], gene_row['gene_id'])).fetchall()

    return jsonify({'coefs': [
        {'tf_name': r['tf_name'], 'grna_name': r['grna_name'], 'active': r['active'],
         'coef': r['coef'], 'signed_pct_rank': r['signed_pct_rank'],
         'has_binding': bool(r['has_binding'])}
        for r in rows
    ]})


@perturbseq_bp.route("/api/module/<mod_name>/enrichment")
def api_module_enrichment(mod_name):
    db = get_db()
    mod_row = db.execute(
        "SELECT module_id FROM module_table WHERE module_name=?", (mod_name,)).fetchone()
    if not mod_row:
        return jsonify([])
    rows = rows_to_dicts(db.execute(
        "SELECT term_id, term_name, source, p_value, term_size, intersection_size "
        "FROM go_module_enrichment WHERE module_id=? AND source = 'GO:BP'ORDER BY p_value",
        (mod_row['module_id'],)).fetchall())
    return jsonify(rows)


@perturbseq_bp.route("/api/search")
def api_search():
    q = request.args.get("q", "").strip()
    kind = request.args.get("kind", "auto")
    if not q:
        return jsonify([])
    db = get_db()
    try:
        if kind == "module":
            rows = db.execute(
                "SELECT module_name FROM module_table "
                "WHERE module_name LIKE ? AND source='hotspot_supermodule' AND module_name != 'unassigned' "
                "ORDER BY module_name LIMIT 20",
                (q + "%",)).fetchall()
            return jsonify([r[0] for r in rows])

        results = []

        # 1. Active TFs (in perturbation library)
        tf_rows = db.execute(
            "SELECT gene_name FROM gene_table WHERE gene_name LIKE ? AND in_perturbation_library=1 "
            "ORDER BY gene_name LIMIT 8",
            (q + "%",)).fetchall()
        active_tf_set = {r[0] for r in tf_rows}
        results.extend({"name": r[0], "type": "active_tf"} for r in tf_rows)

        # 2. Supermodules
        mod_rows = db.execute(
            "SELECT module_name FROM module_table "
            "WHERE module_name LIKE ? AND source='hotspot_supermodule' AND module_name != 'unassigned' "
            "ORDER BY module_name LIMIT 5",
            (q + "%",)).fetchall()
        results.extend({"name": r[0], "type": "module"} for r in mod_rows)

        # 2b. Submodules
        sub_rows = db.execute(
            "SELECT module_name FROM module_table WHERE module_name LIKE ? AND source='hotspot_submodule' "
            "ORDER BY module_name LIMIT 5",
            (q + "%",)).fetchall()
        results.extend({"name": r[0], "type": "submodule"} for r in sub_rows)

        # 2d. GO terms (search by go_name; build URL with accession string)
        go_rows = db.execute(
            "SELECT go_accession, go_name FROM go_term_table "
            "WHERE go_name LIKE ? AND is_obsolete=0 ORDER BY LENGTH(go_name) LIMIT 5",
            ('%' + q + '%',)).fetchall()
        for r in go_rows:
            term_url_name = f"{r['go_name']} ({r['go_accession']})"
            results.append({"name": term_url_name, "type": "go_term",
                            "url": PERTURBSEQ_PREFIX + "/go/" + term_url_name})

        # 2c. TC gene clusters from consolidated DB
        q_up = q.upper()
        for r in db.execute(
                "SELECT module_name FROM module_table WHERE source='mfuzz_k7' ORDER BY module_name"
        ).fetchall():
            internal_id = _mfuzz_to_internal(r[0])
            if internal_id in HIDDEN_TC: continue
            dn = _tc_display(internal_id)
            if dn.upper().startswith(q_up):
                results.append({"name": dn, "type": "gc_cluster",
                                "url": PERTURBSEQ_PREFIX + "/module/" + dn})
        # peak cluster data not available in consolidated DB — omitted

        # 3. Lambert TFs not in the perturbation library
        all_lambert = _load_lambert_tfs()
        q_upper = q.upper()
        lambert_matches = sorted(
            name for name in all_lambert
            if name.upper().startswith(q_upper) and name not in active_tf_set
        )[:8]
        results.extend({"name": name, "type": "tf"} for name in lambert_matches)

        # 4. Non-TF genes
        gene_rows = db.execute(
            "SELECT gene_name FROM gene_table WHERE gene_name LIKE ? ORDER BY gene_name LIMIT 20",
            (q + "%",)).fetchall()
        for r in gene_rows:
            if r[0] not in active_tf_set and r[0] not in all_lambert:
                results.append({"name": r[0], "type": "gene"})
            if len(results) >= 20:
                break

        return jsonify(results[:20])
    except Exception:
        pass
    return jsonify([])


# ── TF-Gene link detail page ──────────────────────────────────────────────────

def _classify_evidence_type(mao_link_type: str | None) -> str:
    if mao_link_type in ('HiCAR', 'HiCAR+TSS'):
        return 'multiome_hicar'
    if mao_link_type == 'TSS':
        return 'multiome'
    return 'proximity'


_ALL_LINK_TYPES = {'proximity', 'multiome', 'multiome_hicar'}


def _parse_link_filter():
    """Read max_distance (int, bp) and link_types (set) from request.args.
    No params present => no filtering (matches current behavior)."""
    max_distance = request.args.get('max_distance', type=int)
    raw = request.args.get('link_types', '')
    link_types = set(raw.split(',')) & _ALL_LINK_TYPES if raw else set(_ALL_LINK_TYPES)
    return max_distance, link_types


def _apg_filter_sql(max_distance, link_types):
    """SQL fragment (leading ' AND ...') + params list to filter an already-aliased
    apg.distance_bp / apg.mao_lt pair. Call after building the apg join."""
    conds, params = [], []
    if max_distance is not None:
        conds.append("apg.distance_bp <= ?")
        params.append(max_distance)
    if link_types != _ALL_LINK_TYPES:
        type_conds = []
        if 'proximity' in link_types:      type_conds.append("apg.mao_lt IS NULL")
        if 'multiome' in link_types:       type_conds.append("apg.mao_lt = 'TSS'")
        if 'multiome_hicar' in link_types: type_conds.append("apg.mao_lt IN ('HiCAR','HiCAR+TSS')")
        conds.append("(" + " OR ".join(type_conds) + ")" if type_conds else "0")
    return (" AND " + " AND ".join(conds)) if conds else "", params


def _classify_distance_label(distance_bp: int | None) -> str:
    if distance_bp is None or distance_bp > 10000:
        return 'distal'
    if distance_bp <= 1000:
        return 'at_tss'
    return 'proximal'


def _norm_chr(chrom: str) -> str:
    """Ensure chromosome has 'chr' prefix (raw data uses bare numbers)."""
    if chrom and not chrom.startswith('chr'):
        return 'chr' + chrom
    return chrom


def _format_tss_dist(dist: int | None) -> str | None:
    """Human-readable distance label for ATAC peak → TSS distance values."""
    if dist is None:
        return None
    if dist == 0:
        return "at TSS"
    if dist >= 1000:
        return f"{dist / 1000:.1f} kb"
    return f"{dist:,} bp"


def _query_tf_gene_link(db, tf: str, gene: str) -> list:
    """Return list of regulatory element dicts for a (TF, gene) pair."""
    gene_row = db.execute("SELECT gene_id FROM gene_table WHERE gene_name=? LIMIT 1", (gene,)).fetchone()
    if not gene_row:
        return []
    gene_id = gene_row['gene_id']
    ds_ids = _tf_dataset_ids(db, tf)
    if not ds_ids:
        return []
    ds_ph = ','.join('?' * len(ds_ids))
    rows = db.execute(f"""
        SELECT
            ap.atac_peak_id, ap.chr, ap.chrom_start AS start, ap.chrom_end AS end,
            tp.source  AS peak_source,
            td.dataset_id, td.dataset, td.cell_type,
            apg.distance_bp, apg.mao_lt
        FROM tf_peaks tp
        JOIN tf_dataset_table td ON td.dataset_id = tp.dataset_id
        JOIN atac_tf_overlaps ato ON ato.peak_id = tp.peak_id
        JOIN atac_peak_table  ap  ON ap.atac_peak_id = ato.atac_peak_id
        JOIN (SELECT atac_peak_id, gene_id, distance_bp, NULL AS mao_lt FROM atac_tss_links
              UNION ALL
              SELECT atac_peak_id, gene_id, distance_to_tss, link_type FROM multiome_atac_overlaps) apg
          ON apg.atac_peak_id = ap.atac_peak_id AND apg.gene_id = ?
        WHERE tp.dataset_id IN ({ds_ph})
    """, [gene_id] + ds_ids).fetchall()

    by_peak: dict = {}
    for r in rows:
        pid = r['atac_peak_id']
        if pid not in by_peak:
            chrom = _norm_chr(str(r['chr']))
            by_peak[pid] = {
                'atac_peak_id': pid,
                'chr':   chrom,
                'start': r['start'],
                'end':   r['end'],
                'element_id': f"{chrom}:{r['start']:,}–{r['end']:,}",
                'motif':   None,
                'chipseq': [],
                '_seen_chips': set(),
                '_seen_motifs': set(),
                '_seen_links':  set(),
                'e2g_links': [],
            }
        el = by_peak[pid]

        if r['peak_source'] == 'tobias':
            key = r['dataset']
            if key not in el['_seen_motifs']:
                el['_seen_motifs'].add(key)
                if el['motif'] is None:
                    el['motif'] = {'name': key, 'score': None, 'pvalue': None}
        else:
            key = r['dataset_id']
            if key not in el['_seen_chips']:
                el['_seen_chips'].add(key)
                el['chipseq'].append({
                    'dataset_id': str(r['dataset_id']),
                    'cell_type':  r['cell_type'],
                    'accession':  r['dataset'],
                })

        link_key = _classify_evidence_type(r['mao_lt'])
        if link_key not in el['_seen_links']:
            el['_seen_links'].add(link_key)
            el['e2g_links'].append({
                'evidence_type':  _classify_evidence_type(r['mao_lt']),
                'distance_label': _classify_distance_label(r['distance_bp']),
                'distance_to_tss': r['distance_bp'],
                'correlation':    None,
                'padj':           None,
                'hicar_score':    None,
                'hicar_padj':     None,
            })

    # Strip internal bookkeeping keys, annotate TSS distance, sort closest-first
    elements = []
    for el in by_peak.values():
        el.pop('_seen_chips', None)
        el.pop('_seen_motifs', None)
        el.pop('_seen_links', None)
        dists = [lnk['distance_to_tss'] for lnk in el['e2g_links'] if lnk['distance_to_tss'] is not None]
        tss_dist = min(dists) if dists else None
        el['tss_dist'] = tss_dist
        el['tss_dist_label'] = _format_tss_dist(tss_dist)
        elements.append(el)
    elements.sort(key=lambda e: (e['tss_dist'] if e['tss_dist'] is not None else 10**9))
    return elements


def _query_tf_binding_peaks(db, tf: str, gene: str) -> list:
    gene_row = db.execute("SELECT gene_id FROM gene_table WHERE gene_name=? LIMIT 1", (gene,)).fetchone()
    if not gene_row:
        return []
    gene_id = gene_row['gene_id']
    ds_ids = _tf_dataset_ids(db, tf)
    if not ds_ids:
        return []
    ds_ph = ','.join('?' * len(ds_ids))
    rows = db.execute(f"""
        SELECT DISTINCT
            tp.chr, tp.chrom_start AS start, tp.chrom_end AS end, tp.tf_peak_name,
            tp.source, td.dataset, td.cell_type, td.dataset_id
        FROM tf_peaks tp
        JOIN tf_dataset_table td ON td.dataset_id = tp.dataset_id
        JOIN atac_tf_overlaps ato ON ato.peak_id = tp.peak_id
        JOIN (SELECT atac_peak_id FROM atac_tss_links WHERE gene_id = ?
              UNION SELECT atac_peak_id FROM multiome_atac_overlaps WHERE gene_id = ?) gp
          ON gp.atac_peak_id = ato.atac_peak_id
        WHERE tp.dataset_id IN ({ds_ph})
    """, [gene_id, gene_id] + ds_ids).fetchall()
    return [
        {
            "chr":       _norm_chr(r["chr"]),
            "start":     r["start"],
            "end":       r["end"],
            "name":      r["tf_peak_name"] or r["dataset"],
            "source":    r["source"],
            "cell_type": r["cell_type"],
            "dataset":   r["dataset"],
        }
        for r in rows
    ]


def _query_element_tfs(db, atac_peak_id: int) -> list:
    rows = db.execute("""
        SELECT DISTINCT
            td.tf_gene_name, td.dataset_id, td.cell_type, td.source AS ds_source,
            td.dataset, tp.source AS peak_source,
            tp.chr, tp.chrom_start AS tf_start, tp.chrom_end AS tf_end,
            ato.overlap_bp
        FROM atac_tf_overlaps ato
        JOIN tf_peaks    tp ON ato.peak_id    = tp.peak_id
        JOIN tf_dataset_table td ON tp.dataset_id  = td.dataset_id
        WHERE ato.atac_peak_id = ?
        ORDER BY td.tf_gene_name, td.cell_type
    """, (atac_peak_id,)).fetchall()
    groups: dict = {}
    for r in rows:
        tf = r["tf_gene_name"]
        if tf not in groups:
            groups[tf] = {"tf": tf, "datasets": []}
        groups[tf]["datasets"].append({
            "dataset_id": r["dataset_id"],
            "cell_type":  r["cell_type"],
            "source":     r["peak_source"],
            "dataset":    r["dataset"],
            "overlap_bp": r["overlap_bp"],
            "tf_chr":     _norm_chr(r["chr"]),
            "tf_start":   r["tf_start"],
            "tf_end":     r["tf_end"],
        })
    return list(groups.values())


def _query_element_genes(db, atac_peak_id: int) -> list:
    rows = db.execute("""
        SELECT DISTINCT
            gtt.gene_name, gtt.gene_id,
            apg.distance_bp, apg.mao_lt,
            td.tf_gene_name,
            gts.tss_position AS tss
        FROM (SELECT atac_peak_id, gene_id, distance_bp, NULL AS mao_lt FROM atac_tss_links
              UNION ALL
              SELECT atac_peak_id, gene_id, distance_to_tss, link_type FROM multiome_atac_overlaps) apg
        JOIN gene_table        gtt ON apg.gene_id      = gtt.gene_id
        LEFT JOIN gene_tss_table gts ON gts.gene_id   = apg.gene_id
        LEFT JOIN atac_tf_overlaps ato ON ato.atac_peak_id = apg.atac_peak_id
        LEFT JOIN tf_peaks         tp  ON ato.peak_id      = tp.peak_id
        LEFT JOIN tf_dataset_table td  ON tp.dataset_id    = td.dataset_id
        WHERE apg.atac_peak_id = ?
    """, (atac_peak_id,)).fetchall()

    groups: dict = {}   # gene_name -> entry dict
    for r in rows:
        g = r["gene_name"]
        dist = r["distance_bp"]
        if g not in groups:
            groups[g] = {
                "gene":          g,
                "mediating_tfs": set(),
                "tss_dist":      dist,
                "tss_dist_label": _format_tss_dist(dist),
                "evidence_type":  _classify_evidence_type(r["mao_lt"]),
                "distance_label": _classify_distance_label(dist),
                "gr_chr":        None,
                "gr_start":      None,
                "gr_end":        None,
                "tss":           r["tss"],
                "tss_chr":       None,
            }
        if r["tf_gene_name"]:
            groups[g]["mediating_tfs"].add(r["tf_gene_name"])
    result = [dict(g, mediating_tfs=sorted(t for t in g["mediating_tfs"] if t is not None)) for g in groups.values()]
    result.sort(key=lambda g: (g["tss_dist"] if g["tss_dist"] is not None else 10**9))
    return result


def _aggregate_atac_counts(rows) -> list:
    """Aggregate a list of atac_peak_counts rows into per-timepoint dicts."""
    by_tp: dict = {}
    for r in rows:
        sample = r['sample_name'] if 'sample_name' in r.keys() else r['sample']
        tp = re.sub(r'_\d+$', '', sample)
        by_tp.setdefault(tp, {'counts': [], 'zscores': []})
        if r['normalized_count'] is not None:
            by_tp[tp]['counts'].append(r['normalized_count'])
        if r['zscore'] is not None:
            by_tp[tp]['zscores'].append(r['zscore'])
    result = []
    for tp in _TP_ORDER:
        if tp not in by_tp:
            continue
        counts  = by_tp[tp]['counts']
        zscores = by_tp[tp]['zscores']
        mean_c  = sum(counts)  / len(counts)  if counts  else None
        mean_z  = sum(zscores) / len(zscores) if zscores else None
        result.append({
            'timepoint':   tp,
            'mean':        round(mean_c, 4) if mean_c is not None else None,
            'mean_z':      round(mean_z, 4) if mean_z is not None else None,
            'replicates':  [round(v, 4) for v in counts],
            'z_replicates': [round(v, 4) for v in zscores],
        })
    return result


def _query_atac_counts(db, atac_peak_id: int) -> list:
    """Return per-timepoint mean normalized_count + z-score for a single peak."""
    try:
        rows = db.execute(
            "SELECT sample_name, normalized_count, zscore FROM atac_peak_counts WHERE atac_peak_id=?",
            (atac_peak_id,)
        ).fetchall()
    except Exception:
        return []
    return _aggregate_atac_counts(rows)



def _query_gene_elements(db, gene: str, gene_id: str | None = None) -> list:
    # Use gene_id (indexed) to avoid a full table scan.
    # Callers that already resolved gene_id should pass it to avoid an extra lookup.
    if gene_id is None:
        row = db.execute(
            "SELECT gene_id FROM gene_table WHERE gene_name=? LIMIT 1", (gene,)).fetchone()
        gene_id = row['gene_id'] if row else None
    if gene_id is None:
        return []

    # Step 1: get distinct ATAC peaks linked to this gene via atac_tss_links or
    # multiome_atac_overlaps, taking the best (smallest distance) link per peak.
    rows = db.execute("""
        SELECT DISTINCT
            ap.atac_peak_id, ap.chr, ap.chrom_start AS start, ap.chrom_end AS end,
            apg.distance_bp, apg.mao_lt
        FROM (SELECT atac_peak_id, gene_id, distance_bp, NULL AS mao_lt FROM atac_tss_links
              UNION ALL
              SELECT atac_peak_id, gene_id, distance_to_tss, link_type FROM multiome_atac_overlaps) apg
        JOIN atac_peak_table ap ON ap.atac_peak_id = apg.atac_peak_id
        WHERE apg.gene_id = ?
        ORDER BY ap.chr, ap.chrom_start
    """, (gene_id,)).fetchall()

    if not rows:
        return []

    # Deduplicate: same ATAC peak can appear via both atac_tss_links and multiome_atac_overlaps.
    # Keep the multiome link if present (it carries HiCAR info), otherwise keep TSS link.
    best: dict = {}  # atac_peak_id -> row
    for r in rows:
        pid = r["atac_peak_id"]
        if pid not in best or (r["mao_lt"] is not None and best[pid]["mao_lt"] is None):
            best[pid] = r
    rows = list(best.values())

    peak_ids = [r["atac_peak_id"] for r in rows]
    placeholders = ",".join("?" * len(peak_ids))

    # Step 2: for each distinct element, count TFs via the atac_peak_id index.
    tf_counts = {
        r["atac_peak_id"]: r["n_tfs"]
        for r in db.execute(f"""
            SELECT ato.atac_peak_id, COUNT(DISTINCT td.tf_gene_name) AS n_tfs
            FROM atac_tf_overlaps ato
            JOIN tf_peaks         tp ON ato.peak_id   = tp.peak_id
            JOIN tf_dataset_table td ON tp.dataset_id = td.dataset_id
            WHERE ato.atac_peak_id IN ({placeholders})
            GROUP BY ato.atac_peak_id
        """, peak_ids).fetchall()
    }

    result = []
    for r in rows:
        dist = r["distance_bp"]
        result.append({
            "atac_peak_id":  r["atac_peak_id"],
            "chr":           _norm_chr(r["chr"]),
            "start":         r["start"],
            "end":           r["end"],
            "evidence_type":  _classify_evidence_type(r["mao_lt"]),
            "distance_label": _classify_distance_label(dist),
            "n_tfs":         tf_counts.get(r["atac_peak_id"], 0),
            "tss_dist":      dist,
            "tss_dist_label": _format_tss_dist(dist),
        })
    return result


@perturbseq_bp.route("/api/link-tfs")
def api_link_tfs():
    db = get_db()
    rows = db.execute(
        "SELECT DISTINCT t.tf_gene_name "
        "FROM tf_dataset_table t "
        "JOIN gene_table g ON g.gene_name = t.tf_gene_name "
        "WHERE t.tf_gene_name IS NOT NULL "
        "  AND g.in_perturbation_library = 1 "
        "ORDER BY t.tf_gene_name"
    ).fetchall()
    return jsonify([r[0] for r in rows])


@perturbseq_bp.route("/api/link-genes")
def api_link_genes():
    db = get_db()
    rows = db.execute(
        "SELECT DISTINCT gt.gene_name "
        "FROM (SELECT gene_id FROM atac_tss_links "
        "      UNION SELECT gene_id FROM multiome_atac_overlaps) apg "
        "JOIN gene_table gt ON apg.gene_id = gt.gene_id "
        "ORDER BY gt.gene_name"
    ).fetchall()
    return jsonify([r[0] for r in rows])


@perturbseq_bp.route("/api/link-modules")
def api_link_modules():
    db = get_db()
    rows = db.execute(
        "SELECT module_name FROM module_table "
        "WHERE source='hotspot_supermodule' AND module_name != 'unassigned' ORDER BY module_name"
    ).fetchall()
    return jsonify([r[0] for r in rows])


@perturbseq_bp.route("/tf-module/<tf>/<module_name>")
def tf_module_link_page(tf, module_name):
    db = get_db()

    tf_row = resolve_gene(db, tf)
    if not tf_row:
        abort(404)
    tf = tf_row['gene_name']

    mod_row = db.execute(
        "SELECT module_id, module_name, size FROM module_table "
        "WHERE module_name=? AND source='hotspot_supermodule'",
        (module_name,)).fetchone()
    if not mod_row:
        abort(404)
    module_id = mod_row['module_id']

    pert_row = db.execute(
        "SELECT mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj, direction "
        "FROM gsea_tf_table "
        "WHERE gene_name=? AND module=? AND module_collection='hotspot_supermodule' "
        "AND gene_set_collection='DE-hotspot_modules'",
        (tf, module_name)).fetchone()
    pert = dict(pert_row) if pert_row else None

    grna_nes = rows_to_dicts(db.execute(
        "SELECT g.grna_id, g.NES, g.padj, gt.gene_name "
        "FROM gsea_grna_table g JOIN grna_table gt ON g.grna_id=gt.grna_id "
        "WHERE gt.gene_name=? AND g.module=? AND g.module_collection='hotspot_supermodule' "
        "ORDER BY g.NES DESC",
        (tf, module_name)))

    binding_datasets = rows_to_dicts(db.execute("""
        SELECT td.source, td.cell_type, td.dataset, te.odds_ratio AS odds_ratio,
               te.padj_fisher, te.n_overlap, te.n_tf_targets
        FROM tf_module_enrichment te
        JOIN tf_dataset_table td ON td.tf_gene_name = te.tf_gene_name
        WHERE te.tf_gene_name=? AND te.module_id=?
          AND te.gene_set_collection='hotspot_supermodule'
        ORDER BY te.padj_fisher
    """, (tf, module_id)))
    best_binding = binding_datasets[0] if binding_datasets else None

    sub_pert = rows_to_dicts(db.execute(
        "SELECT module AS submodule, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj, direction "
        "FROM gsea_tf_table "
        "WHERE gene_name=? AND module LIKE ? AND module_collection='hotspot_submodule' "
        "AND gene_set_collection='DE-hotspot_modules' "
        "ORDER BY ABS(mean_NES) DESC",
        (tf, module_name + '.%')))

    tf_expr    = _gene_expr_by_timepoint(db, tf_row['gene_id'])
    module_expr = _module_expr_by_timepoint(db, module_id)

    go_terms = rows_to_dicts(db.execute(
        "SELECT term_name, p_value, source FROM go_module_enrichment "
        "WHERE module_id=? AND source IN ('GO:BP','GO:MF','GO:CC') "
        "ORDER BY p_value LIMIT 8",
        (module_id,)))

    top_genes = [r[0] for r in db.execute(
        "SELECT gt.gene_name FROM gene_module_table gm "
        "JOIN gene_table gt ON gm.gene_id = gt.gene_id "
        "WHERE gm.module_id=? AND gt.gene_name NOT LIKE 'ENSG%' "
        "ORDER BY gt.gene_name LIMIT 40",
        (module_id,)).fetchall()]

    has_pert = pert is not None
    has_bind = bool(binding_datasets)
    if not has_pert and not has_bind:
        abort(404)
    evidence = 'both' if (has_pert and has_bind) else ('perturbation' if has_pert else 'binding')

    module_color = MODULE_COLORS.get(module_name, '#4a56b0')
    tfs, modules = _nav_lists(db)

    return render_template(
        "perturbseq/tf_module_link.html",
        tf=tf,
        module=module_name,
        module_size=mod_row['size'],
        module_color=module_color,
        pert=pert,
        grna_nes=grna_nes,
        binding_datasets=binding_datasets,
        best_binding=best_binding,
        sub_pert=sub_pert,
        tf_expr=tf_expr,
        module_expr=module_expr,
        go_terms=go_terms,
        top_genes=top_genes,
        evidence=evidence,
        tfs=tfs,
        modules=modules,
    )


def _get_perturbation_evidence(db, tf: str, gene: str) -> dict:
    """Return perturbation evidence for a (TF, gene) pair.

    Returns a dict:
      gene_in_de: bool — whether the gene appears in the DE regression data at all
      hits:       list — {grna, direction, coef} for gRNAs that place gene in top/bottom 5%
                  direction is 'down' (bottom 5%) or 'up' (top 5%)
    """
    rows = db.execute("""
        SELECT gt.grna_id, gt.grna_name, d.coef
        FROM de_results d
        JOIN grna_table gt ON d.grna_id = gt.grna_id
        WHERE gt.gene_name = ?
        ORDER BY gt.grna_id
    """, (tf,)).fetchall()
    if not rows:
        return {'gene_in_de': False, 'hits': []}

    grna_coefs: dict = {}
    grna_names: dict = {}
    for grna_id, grna_name, coef in rows:
        grna_coefs.setdefault(grna_id, []).append(coef)
        grna_names[grna_id] = grna_name

    target_rows = db.execute("""
        SELECT d.grna_id, d.coef
        FROM de_results d
        JOIN grna_table gt ON d.grna_id = gt.grna_id
        JOIN gene_table g ON d.gene_id = g.gene_id
        WHERE gt.gene_name = ? AND g.gene_name = ?
    """, (tf, gene)).fetchall()
    if not target_rows:
        return {'gene_in_de': False, 'hits': []}

    target_map = {r[0]: r[1] for r in target_rows}
    hits = []
    for grna_id, coefs in grna_coefs.items():
        target_coef = target_map.get(grna_id)
        if target_coef is None:
            continue
        arr = np.array(coefs, dtype=float)
        lo5 = float(np.percentile(arr, 5))
        hi5 = float(np.percentile(arr, 95))
        if target_coef <= lo5:
            hits.append({'grna': grna_names[grna_id], 'direction': 'down', 'coef': round(target_coef, 4)})
        elif target_coef >= hi5:
            hits.append({'grna': grna_names[grna_id], 'direction': 'up',   'coef': round(target_coef, 4)})
    return {'gene_in_de': True, 'hits': hits}


def _coef_kde_data_db(db, tf: str, gene: str, n_pts: int = 100):
    """Return KDE density curves for all gRNAs targeting tf, using DE coefficients from DB.

    Returns None if no gRNAs found for this TF.
    """
    rows = db.execute("""
        SELECT gt.grna_id, gt.grna_name, d.coef
        FROM de_results d
        JOIN grna_table gt ON d.grna_id = gt.grna_id
        WHERE gt.gene_name = ?
        ORDER BY gt.grna_id
    """, (tf,)).fetchall()

    if not rows:
        return None

    grna_coefs = defaultdict(list)
    grna_names = {}
    for grna_id, grna_name, coef in rows:
        grna_coefs[grna_id].append(coef)
        grna_names[grna_id] = grna_name

    target_rows = db.execute("""
        SELECT d.grna_id, d.coef
        FROM de_results d
        JOIN grna_table gt ON d.grna_id = gt.grna_id
        JOIN gene_table g ON d.gene_id = g.gene_id
        WHERE gt.gene_name = ? AND g.gene_name = ?
    """, (tf, gene)).fetchall()
    target_map = {r[0]: r[1] for r in target_rows}

    all_vals = [c for coefs in grna_coefs.values() for c in coefs]
    x_min = float(min(all_vals))
    x_max = float(max(all_vals))
    rng = (x_max - x_min) or 1.0
    pad = 0.02 * rng
    xs = np.linspace(x_min - pad, x_max + pad, n_pts)

    result_grnas = []
    for grna_id, coefs in grna_coefs.items():
        coefs_arr = np.array(coefs, dtype=float)
        try:
            ys = gaussian_kde(coefs_arr, bw_method='scott')(xs).tolist()
        except Exception:
            ys = [0.0] * n_pts
        target_coef = target_map.get(grna_id)
        lo5 = float(np.percentile(coefs_arr, 5))
        hi5 = float(np.percentile(coefs_arr, 95))
        result_grnas.append({
            'grna': grna_names[grna_id],
            'y': [round(v, 6) for v in ys],
            'target_coef': round(target_coef, 6) if target_coef is not None else None,
            'lo5': round(lo5, 5),
            'hi5': round(hi5, 5),
        })

    return {'x': [round(v, 5) for v in xs.tolist()], 'gRNAs': result_grnas}


@perturbseq_bp.route("/api/link/<tf>/<gene>/coefs")
def api_link_coefs(tf, gene):
    db = get_db()
    try:
        de = _coef_kde_data_db(db, tf, gene)
    except Exception:
        return jsonify({"tf": tf, "gene": gene, "de": None,
                        "error": "coefficient data unavailable"}), 503
    return jsonify({"tf": tf, "gene": gene, "de": de})


@perturbseq_bp.route("/api/link/<tf>/<gene>/tooltip")
def api_link_tooltip(tf, gene):
    db = get_db()
    tf_row   = resolve_gene(db, tf)
    gene_row = resolve_gene(db, gene)
    tf_expr = [{'timepoint': r['timepoint'], 'mean_tpm': r['mean_tpm']}
               for r in _gene_expr_by_timepoint(db, tf_row['gene_id'])] if tf_row else []
    gene_expr = [{'timepoint': r['timepoint'], 'mean_tpm': r['mean_tpm']}
                 for r in _gene_expr_by_timepoint(db, gene_row['gene_id'])] if gene_row else []
    return jsonify({'tf': tf, 'gene': gene, 'tf_expr': tf_expr, 'gene_expr': gene_expr})


@perturbseq_bp.route("/api/atac-counts/<int:atac_peak_id>")
def atac_counts_api(atac_peak_id):
    """JSON: per-timepoint accessibility profile for a single ATAC peak."""
    db = get_db()
    return jsonify(_query_atac_counts(db, atac_peak_id))


@perturbseq_bp.route("/bw/<path:filename>")
def serve_bw(filename):
    return send_from_directory(BW_DIR, filename)


@perturbseq_bp.route("/bw-extra/<path:filename>")
def serve_bw_extra(filename):
    return send_from_directory(BW_EXTRA_DIR, filename)


@perturbseq_bp.route("/bw-rna/<path:filename>")
def serve_bw_rna(filename):
    return send_from_directory(BW_RNA_DIR, filename)


@perturbseq_bp.route("/gtf/<path:filename>")
def serve_gtf(filename):
    response = send_from_directory(GTF_DIR, filename)
    response.headers.pop('Content-Encoding', None)
    return response


def _build_element_description(peak: dict, genes: list, tfs: list) -> str:
    size_bp = peak["end"] - peak["start"]
    size_str = f"{size_bp / 1000:.1f} kb" if size_bp >= 1000 else f"{size_bp:,} bp"

    n_genes = len(genes)
    if n_genes == 0:
        gene_part = "has no linked target genes"
    elif n_genes == 1:
        gene_part = f"is linked to {genes[0]['gene']}"
    elif n_genes == 2:
        gene_part = f"is linked to {genes[0]['gene']} and {genes[1]['gene']}"
    elif n_genes == 3:
        names = [g["gene"] for g in genes]
        gene_part = f"is linked to {names[0]}, {names[1]}, and {names[2]}"
    else:
        gene_part = f"is linked to {n_genes} target genes"

    n_tfs = len(tfs)
    if n_tfs == 0:
        tf_part = "shows no TF binding evidence"
    elif n_tfs == 1:
        tf_part = f"is bound by {tfs[0]['tf']}"
    elif n_tfs == 2:
        tf_part = f"is bound by {tfs[0]['tf']} and {tfs[1]['tf']}"
    elif n_tfs <= 5:
        names = [t["tf"] for t in tfs]
        tf_part = f"is bound by {', '.join(names[:-1])}, and {names[-1]}"
    else:
        tf_part = f"is bound by {n_tfs} transcription factors"

    return f"This {size_str} open chromatin element {gene_part} and {tf_part}."


def _query_element_gene_link(db, atac_peak_id: int, gene: str) -> dict | None:
    """Return TF-mediated link details for a specific (element, gene) pair."""
    peak_row = db.execute(
        'SELECT atac_peak_id, chr, chrom_start AS start, chrom_end AS end FROM atac_peak_table WHERE atac_peak_id=?',
        (atac_peak_id,)
    ).fetchone()
    if peak_row is None:
        return None
    peak = dict(peak_row)
    peak["chr"] = _norm_chr(peak["chr"])

    gene_row = db.execute("SELECT gene_id FROM gene_table WHERE gene_name=? LIMIT 1", (gene,)).fetchone()
    gene_id = gene_row['gene_id'] if gene_row else None
    if not gene_id:
        return None

    # Determine evidence_type/distance from atac_tss_links / multiome_atac_overlaps
    link_row = db.execute("""
        SELECT distance_bp, NULL AS mao_lt FROM atac_tss_links
        WHERE atac_peak_id = ? AND gene_id = ?
        UNION ALL
        SELECT distance_to_tss, link_type FROM multiome_atac_overlaps
        WHERE atac_peak_id = ? AND gene_id = ?
        ORDER BY distance_bp LIMIT 1
    """, (atac_peak_id, gene_id, atac_peak_id, gene_id)).fetchone()

    if not link_row:
        return None

    dist = link_row["distance_bp"]
    evidence_type  = _classify_evidence_type(link_row["mao_lt"])
    distance_label = _classify_distance_label(dist)
    dist_label     = _format_tss_dist(dist)

    rows = db.execute("""
        SELECT DISTINCT
            td.tf_gene_name,
            td.dataset_id, td.cell_type, td.dataset,
            tp.source AS peak_source,
            ato.overlap_bp
        FROM atac_tf_overlaps ato
        JOIN tf_peaks         tp ON ato.peak_id    = tp.peak_id
        JOIN tf_dataset_table td ON tp.dataset_id  = td.dataset_id
        WHERE ato.atac_peak_id = ?
        ORDER BY td.tf_gene_name, td.cell_type
    """, (atac_peak_id,)).fetchall()

    tfs: dict = {}
    for r in rows:
        tf = r["tf_gene_name"]
        if tf not in tfs:
            tfs[tf] = {"tf": tf, "datasets": [], "_seen": set()}
        key = (r["dataset_id"], r["peak_source"])
        if key not in tfs[tf]["_seen"]:
            tfs[tf]["_seen"].add(key)
            tfs[tf]["datasets"].append({
                "dataset_id": r["dataset_id"],
                "cell_type":  r["cell_type"],
                "source":     r["peak_source"],
                "dataset":    r["dataset"],
                "overlap_bp": r["overlap_bp"],
            })
    for tf_data in tfs.values():
        tf_data.pop("_seen")

    return {
        "peak":         peak,
        "gene":         gene,
        "evidence_type":  evidence_type,
        "distance_label": distance_label,
        "dist_label":   dist_label,
        "mediating_tfs": sorted(tfs.values(), key=lambda t: t["tf"]),
    }


@perturbseq_bp.route("/element/<int:atac_peak_id>/link/<gene>")
def element_gene_link_page(atac_peak_id, gene):
    db   = get_db()
    data = _query_element_gene_link(db, atac_peak_id, gene)
    if data is None:
        abort(404)
    atac_counts = _query_atac_counts(db, atac_peak_id)
    gene_row = resolve_gene(db, gene)
    gene_expr = [{'timepoint': r['timepoint'], 'mean_tpm': r['mean_tpm']}
                 for r in _gene_expr_by_timepoint(db, gene_row['gene_id'])] if gene_row else []
    return render_template(
        "perturbseq/element_gene_link.html",
        peak=data["peak"],
        gene=gene,
        evidence_type=data["evidence_type"],
        distance_label=data["distance_label"],
        dist_label=data["dist_label"],
        mediating_tfs=data["mediating_tfs"],
        atac_counts=atac_counts,
        gene_expr=gene_expr,
    )


@perturbseq_bp.route("/element/<int:atac_peak_id>")
def element_page(atac_peak_id):
    db = get_db()
    row = db.execute(
        'SELECT atac_peak_id, chr, chrom_start AS start, chrom_end AS end, atac_peak_name FROM atac_peak_table WHERE atac_peak_id=?',
        (atac_peak_id,)
    ).fetchone()
    if row is None:
        abort(404)
    peak = dict(row)
    peak["chr"] = _norm_chr(peak["chr"])

    tfs         = _query_element_tfs(db, atac_peak_id)
    genes       = _query_element_genes(db, atac_peak_id)
    atac_counts = _query_atac_counts(db, atac_peak_id)

    # Build flat TF peak list for IGV (all datasets with coordinates)
    tf_peaks_igv = []
    for tf_group in tfs:
        for ds in tf_group["datasets"]:
            if ds["tf_chr"] and ds["tf_start"] is not None:
                tf_peaks_igv.append({
                    "chr":    ds["tf_chr"],
                    "start":  ds["tf_start"],
                    "end":    ds["tf_end"],
                    "name":   f"{tf_group['tf']} · {ds['cell_type']}",
                    "source": ds["source"],
                })

    description = _build_element_description(peak, genes, tfs)

    # Nearby ATAC peaks for the "other elements" IGV track (±500 kb window)
    mid = (peak["start"] + peak["end"]) // 2
    nearby_rows = db.execute(
        """SELECT chr, chrom_start AS start, chrom_end AS end, atac_peak_name
           FROM atac_peak_table
           WHERE chr = ? AND chrom_start >= ? AND chrom_end <= ? AND atac_peak_id != ?
           ORDER BY chrom_start""",
        (str(peak["chr"]).lstrip("chr"), mid - 500000, mid + 500000, atac_peak_id)
    ).fetchall()
    nearby_peaks = [
        {"chr": peak["chr"], "start": r["start"], "end": r["end"], "name": r["atac_peak_name"] or ""}
        for r in nearby_rows
    ]

    return render_template(
        "perturbseq/element.html",
        peak=peak,
        tfs=tfs,
        genes=genes,
        tf_peaks_igv=tf_peaks_igv,
        atac_counts=atac_counts,
        nearby_peaks=nearby_peaks,
        description=description,
    )


@perturbseq_bp.route("/link/<tf>/<gene>")
def tf_gene_link_page(tf, gene):
    db = get_db()

    elements = _query_tf_gene_link(db, tf, gene)
    if elements:
        peak_ids = [el['atac_peak_id'] for el in elements]
        ph = ','.join('?' * len(peak_ids))
        count_rows = db.execute(
            f"SELECT atac_peak_id, sample_name, normalized_count, zscore "
            f"FROM atac_peak_counts WHERE atac_peak_id IN ({ph})",
            peak_ids).fetchall()
        counts_by_peak: dict = {}
        for r in count_rows:
            counts_by_peak.setdefault(r['atac_peak_id'], []).append(r)
        for el in elements:
            el['atac_counts'] = _aggregate_atac_counts(counts_by_peak.get(el['atac_peak_id'], []))
    tf_peaks_data = _query_tf_binding_peaks(db, tf, gene)

    tf_row = resolve_gene(db, tf)
    gene_row_data = resolve_gene(db, gene)
    tf_expr = [{'timepoint': r['timepoint'], 'mean_tpm': r['mean_tpm']}
               for r in _gene_expr_by_timepoint(db, tf_row['gene_id'])] if tf_row else []
    gene_expr = [{'timepoint': r['timepoint'], 'mean_tpm': r['mean_tpm']}
                 for r in _gene_expr_by_timepoint(db, gene_row_data['gene_id'])] if gene_row_data else []

    def _get_summary(row):
        if not row:
            return None
        r = db.execute("SELECT Summary FROM gene_table WHERE gene_id=?", (row['gene_id'],)).fetchone()
        return r['Summary'] if r else None

    gene_tss = None
    if gene_row_data:
        _tss = db.execute(
            "SELECT tss_position AS tss, chr FROM gene_tss_table WHERE gene_id=? ORDER BY tss_id LIMIT 1",
            (gene_row_data['gene_id'],)).fetchone()
        if _tss:
            gene_tss = {'tss': _tss['tss'], 'chr': _norm_chr(str(_tss['chr']))}

    has_binding_evidence      = bool(elements)
    _pert                     = _get_perturbation_evidence(db, tf, gene)
    perturbation_evidence     = _pert['hits']
    gene_in_de                = _pert['gene_in_de']
    has_perturbation_evidence = bool(perturbation_evidence)

    return render_template(
        "perturbseq/tf_gene_link.html",
        tf=tf,
        gene=gene,
        tf_desc=_get_summary(tf_row),
        gene_desc=_get_summary(gene_row_data),
        elements=elements,
        tf_peaks_data=tf_peaks_data,
        tf_expr=tf_expr,
        gene_expr=gene_expr,
        gene_tss=gene_tss,
        has_binding_evidence=has_binding_evidence,
        has_perturbation_evidence=has_perturbation_evidence,
        perturbation_evidence=perturbation_evidence,
        gene_in_de=gene_in_de,
    )


@perturbseq_bp.route("/dataset/<int:dataset_id>")
def dataset_page(dataset_id):
    db = get_db()

    ds = db.execute(
        "SELECT * FROM tf_dataset_table WHERE dataset_id=?", (dataset_id,)
    ).fetchone()
    if not ds:
        abort(404)

    tf_summary_row = db.execute(
        "SELECT Summary FROM gene_table WHERE gene_name=?", (ds["tf_gene_name"],)
    ).fetchone()
    tf_desc = {'name': ds["tf_gene_name"],
               'summary': tf_summary_row['Summary'] if tf_summary_row else None}

    stats = {
        "peak_count": db.execute(
            "SELECT COUNT(*) FROM tf_peaks WHERE dataset_id=?", (dataset_id,)
        ).fetchone()[0],
        "elements_overlapped": db.execute(
            """SELECT COUNT(DISTINCT ato.atac_peak_id)
               FROM atac_tf_overlaps ato
               JOIN tf_peaks tp ON ato.peak_id = tp.peak_id
               WHERE tp.dataset_id=?""", (dataset_id,)
        ).fetchone()[0],
        "target_genes": db.execute(
            "SELECT COUNT(DISTINCT gene_id) FROM tf_dataset_gene_table WHERE dataset_id=?",
            (dataset_id,)
        ).fetchone()[0],
    }

    gene_rows = rows_to_dicts(db.execute(
        "SELECT gt.gene_name, dg.match_type FROM tf_dataset_gene_table dg "
        "JOIN gene_table gt ON dg.gene_id = gt.gene_id "
        "WHERE dg.dataset_id=? ORDER BY gt.gene_name",
        (dataset_id,)))

    element_rows = rows_to_dicts(db.execute("""
        SELECT ap.atac_peak_id, ap.chr, ap.chrom_start AS start, ap.chrom_end AS end,
               MAX(ato.overlap_bp) AS overlap_bp,
               COUNT(DISTINCT apg.gene_id) AS gene_count
        FROM atac_tf_overlaps ato
        JOIN tf_peaks        tp ON ato.peak_id      = tp.peak_id
        JOIN atac_peak_table ap ON ato.atac_peak_id = ap.atac_peak_id
        LEFT JOIN (SELECT atac_peak_id, gene_id FROM atac_tss_links
                   UNION ALL
                   SELECT atac_peak_id, gene_id FROM multiome_atac_overlaps) apg
          ON apg.atac_peak_id = ato.atac_peak_id
        WHERE tp.dataset_id=?
        GROUP BY ap.atac_peak_id
        ORDER BY ap.chr, ap.chrom_start
    """, (dataset_id,)))
    for r in element_rows:
        r["chr"] = _norm_chr(str(r["chr"]))

    related = rows_to_dicts(db.execute("""
        SELECT dataset_id, dataset, cell_type, cell_type_group, source
        FROM tf_dataset_table
        WHERE tf_gene_name=? AND dataset_id!=?
        ORDER BY cell_type, dataset
    """, (ds["tf_gene_name"], dataset_id)))

    return render_template(
        "perturbseq/dataset.html",
        ds=ds,
        tf_desc=tf_desc,
        stats=stats,
        gene_rows=gene_rows,
        element_rows=element_rows,
        related=related,
    )


@perturbseq_bp.route("/db-schema")
def db_schema():
    db = get_db()
    tables = db.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    ).fetchall()
    schema = {}
    fks = {}
    for t in tables:
        name = t["name"]
        cols = db.execute(f"PRAGMA table_info({name})").fetchall()
        schema[name] = [
            {"cid": c["cid"], "name": c["name"], "type": c["type"], "pk": bool(c["pk"])}
            for c in cols
        ]
        fk_rows = db.execute(f"PRAGMA foreign_key_list({name})").fetchall()
        fks[name] = [
            {"from_col": f["from"], "to_table": f["table"], "to_col": f["to"]}
            for f in fk_rows
        ]
    return render_template("perturbseq/db_schema.html", schema=schema, fks=fks, db_path=DB_PATH)
