import csv
import json
import math
import re
import sqlite3
from collections import defaultdict
from pathlib import Path
from flask import Blueprint, render_template, jsonify, request, g, redirect, url_for, abort, send_from_directory

perturbseq_bp = Blueprint('perturbseq', __name__, url_prefix='/perturbseq')

DB_PATH = str(Path(__file__).resolve().parent.parent.parent / "data" / "perturbseq.db")
BW_DIR       = "/home/torred1/pipelines/tf-perturbseq/data/bulk-atacseq/s08-bigwig.dir"
BW_EXTRA_DIR = "/data1/huangfud/torred1/sandbox/sandbox008-hotspot_network/data/08-2026_05_27_data/bw"
BW_RNA_DIR   = "/data1/huangfud/torred1/sandbox/sandbox008-hotspot_network/data/09-2026_05_28_data/bw"
GTF_DIR      = "/data1/huangfud/torred1/sandbox/sandbox008-hotspot_network/data/08-2026_05_27_data/gtf"
TC_DB_PATH = str(Path(__file__).resolve().parent / "timecourse.db")

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

def _top_genes_by_pub(genes: list, n: int = 8) -> list:
    pub_counts = _load_gene_pub_counts()
    return sorted(genes, key=lambda g: pub_counts.get(g, 0), reverse=True)[:n]


def _compute_mean_zscore_profile(db, genes: list, timepoints: list) -> list:
    """Z-score each gene's trajectory (row-normalise), return mean ± SEM per timepoint."""
    if not genes or not timepoints:
        return []
    ph = ','.join('?' * len(genes))
    rows = db.execute(
        f"SELECT gene, timepoint, mean_tpm FROM gene_expression WHERE gene IN ({ph})",
        genes).fetchall()
    gene_tpm: dict = defaultdict(dict)
    for gene, tp, tpm in rows:
        gene_tpm[gene][tp] = tpm
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
    ph = ','.join('?' * len(genes))
    rows = db.execute(
        f"SELECT timepoint, mean_tpm FROM gene_expression WHERE gene IN ({ph})",
        genes).fetchall()
    tp_vals: dict = defaultdict(list)
    for tp, tpm in rows:
        if tp in timepoints and tpm is not None:
            tp_vals[tp].append(tpm)
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


def _tc_display(cid: str) -> str:
    if cid.startswith('gene_cluster_'):
        return 'GC' + cid[len('gene_cluster_'):]
    if cid.startswith('peak_cluster_'):
        return 'PC' + cid[len('peak_cluster_'):]
    return cid

_expressed_tfs: frozenset | None = None
_lambert_tfs: frozenset | None = None
_perturbed_tfs: frozenset | None = None
_gene_name_map: dict | None = None
_gene_pub_counts: dict | None = None
_submodule_descriptions: dict | None = None

def _load_lambert_tfs() -> frozenset:
    global _lambert_tfs
    if _lambert_tfs is None:
        try:
            path = Path(__file__).resolve().parent.parent.parent / "data" / "tables" / "lambert-TF_table.csv"
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
            path = Path(__file__).resolve().parent.parent.parent / "data" / "tables" / "gene_pub_counts.tsv"
            counts = {}
            with open(path) as f:
                for line in f:
                    parts = line.rstrip('\n').split('\t', 1)
                    if len(parts) == 2:
                        counts[parts[0]] = int(parts[1])
            _gene_pub_counts = counts
        except Exception:
            _gene_pub_counts = {}
    return _gene_pub_counts


def _load_gene_name_map() -> dict:
    global _gene_name_map
    if _gene_name_map is None:
        try:
            path = Path(__file__).resolve().parent / "gene_name_map.tsv"
            with open(path) as f:
                _gene_name_map = {}
                for line in f:
                    parts = line.rstrip('\n').split('\t', 1)
                    if len(parts) == 2:
                        _gene_name_map[parts[0]] = parts[1]
        except Exception:
            _gene_name_map = {}
    return _gene_name_map


def _load_submodule_descriptions() -> dict:
    """Load LLM-generated submodule descriptions from TSV."""
    global _submodule_descriptions
    if _submodule_descriptions is None:
        try:
            path = Path(__file__).resolve().parent.parent.parent / "data" / "06-llm_submodule_descriptions" / "submodule_descriptions-v3.tsv"
            descs = {}
            with open(path) as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    descs[row['submodule']] = row['description']
            _submodule_descriptions = descs
        except Exception:
            _submodule_descriptions = {}
    return _submodule_descriptions


def get_tc_db():
    if "tc_db" not in g:
        g.tc_db = sqlite3.connect(TC_DB_PATH)
        g.tc_db.row_factory = sqlite3.Row
    return g.tc_db


def _load_expressed_tfs() -> frozenset:
    """Load expressed TF set from DB once per worker process."""
    global _expressed_tfs
    if _expressed_tfs is None:
        try:
            with sqlite3.connect(DB_PATH) as conn:
                rows = conn.execute("SELECT tf FROM expressed_tfs").fetchall()
            _expressed_tfs = frozenset(r[0] for r in rows)
        except Exception:
            _expressed_tfs = frozenset()  # fail open if table not yet populated
    return _expressed_tfs


def _load_perturbed_tfs() -> frozenset:
    """TFs with at least one significant perturbation association (appear in perturbation_edges)."""
    global _perturbed_tfs
    if _perturbed_tfs is None:
        try:
            with sqlite3.connect(DB_PATH) as conn:
                rows = conn.execute("SELECT DISTINCT tf FROM perturbation_edges").fetchall()
            _perturbed_tfs = frozenset(r[0] for r in rows)
        except Exception:
            _perturbed_tfs = frozenset()
    return _perturbed_tfs


def get_db():
    if "perturbseq_db" not in g:
        g.perturbseq_db = sqlite3.connect(DB_PATH)
        g.perturbseq_db.row_factory = sqlite3.Row
    return g.perturbseq_db


@perturbseq_bp.context_processor
def inject_globals():
    return {'module_colors': MODULE_COLORS}


@perturbseq_bp.teardown_request
def close_db(exc):
    db = g.pop("perturbseq_db", None)
    if db is not None:
        db.close()
    tc = g.pop("tc_db", None)
    if tc is not None:
        tc.close()


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


@perturbseq_bp.route("/")
def index():
    db = get_db()
    tfs = [r[0] for r in db.execute(
        "SELECT DISTINCT tf FROM perturbation_edges UNION SELECT DISTINCT tf FROM binding_edges ORDER BY 1")]
    modules = [r[0] for r in db.execute(
        "SELECT DISTINCT module FROM perturbation_edges WHERE module_type='supermodule' "
        "UNION SELECT DISTINCT module FROM binding_edges WHERE module_type='supermodule' ORDER BY 1")]
    return render_template("perturbseq/index.html", tfs=tfs, modules=modules)


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


def _tp_order():
    return ("ORDER BY CASE timepoint "
            "WHEN 'ES_0h' THEN 1 WHEN 'DE_12h' THEN 2 WHEN 'DE_24h' THEN 3 "
            "WHEN 'DE_36h' THEN 4 WHEN 'DE_48h' THEN 5 WHEN 'DE_60h' THEN 6 "
            "WHEN 'DE_72h' THEN 7 ELSE 8 END")


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
            'evidence': 'both' if has_bind else 'perturbation',
        })
    for k, b in bind_map.items():
        if k not in seen:
            rows.append({
                id_key: k, 'direction': '—', 'mean_NES': None,
                'n_sig_gRNA': None, 'padj': None, 'binding': '✓',
                'odds_ratio': b['odds_ratio'], 'evidence': 'binding',
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


def _nav_lists(db):
    tfs = [r[0] for r in db.execute("SELECT DISTINCT tf FROM tf_descriptions ORDER BY 1")]
    modules = [r[0] for r in db.execute(
        "SELECT DISTINCT module FROM module_descriptions WHERE module NOT LIKE 'GC%' ORDER BY 1")]
    return tfs, modules


def _super_page(mod_name):
    db = get_db()
    exists = db.execute("SELECT 1 FROM module_genes WHERE module=? LIMIT 1", (mod_name,)).fetchone()
    if not exists:
        return render_template("perturbseq/404.html", message=f"Module not found: {mod_name}"), 404
    desc_row = db.execute("SELECT * FROM module_descriptions WHERE module=?", (mod_name,)).fetchone()
    desc = None
    if desc_row:
        desc = dict(desc_row)
        try:
            desc["highlight_terms"] = json.loads(desc["highlight_terms"])
        except (json.JSONDecodeError, TypeError):
            desc["highlight_terms"] = []
        try:
            desc["top_tfs"] = json.loads(desc["top_tfs"])
        except (json.JSONDecodeError, TypeError):
            desc["top_tfs"] = []
    expr = rows_to_dicts(db.execute(
        f"SELECT timepoint, mean_tpm, sd_tpm, n_genes FROM module_expression WHERE module=? {_tp_order()}",
        (mod_name,)))
    genes_list = [r[0] for r in db.execute(
        "SELECT gene FROM module_genes WHERE module=?", (mod_name,))]
    if desc:
        desc["highlight_genes"] = _top_genes_by_pub(genes_list, n=5)
    timepoints = [r['timepoint'] for r in expr]
    expr_z = _compute_mean_zscore_profile(db, genes_list, timepoints)
    submodules = rows_to_dicts(db.execute(
        "SELECT id, n_genes, color, within_mean_z FROM submodule_nodes WHERE supermodule=? ORDER BY n_genes DESC",
        (mod_name,)))
    if submodules:
        pub_counts = _load_gene_pub_counts()
        sub_ids = [s['id'] for s in submodules]
        ph = ','.join('?' * len(sub_ids))
        sub_gene_rows = db.execute(
            f"SELECT submodule, gene FROM submodule_genes WHERE submodule IN ({ph})",
            sub_ids).fetchall()
        sub_genes_map: dict = defaultdict(list)
        for submod, gene in sub_gene_rows:
            sub_genes_map[submod].append(gene)
        for s in submodules:
            genes = sub_genes_map.get(s['id'], [])
            s['top_genes'] = sorted(genes, key=lambda g: pub_counts.get(g, 0), reverse=True)[:5]
    tfs, modules = _nav_lists(db)
    mod_color = MODULE_COLORS.get(mod_name, '#4a56b0')
    return render_template("perturbseq/module.html",
        module_type='super', name=mod_name, color=mod_color, card_color=mod_color,
        desc=desc, expr=expr, expr_z=expr_z, profile=None, parent=None, node=None,
        submodules=submodules, neighbors=[], genes=genes_list,
        notable_genes=[],
        top_tfs=[], assoc=[], assoc_label='', enrichment=[],
        tfs=tfs, modules=modules)


def _sub_page(name):
    db = get_db()
    node_row = db.execute("SELECT * FROM submodule_nodes WHERE id=?", (name,)).fetchone()
    if not node_row:
        return render_template("perturbseq/404.html", message=f"Submodule not found: {name}"), 404
    node = dict(node_row)
    genes = [r[0] for r in db.execute(
        "SELECT gene FROM submodule_genes WHERE submodule=? ORDER BY gene", (name,))]
    pub_counts = _load_gene_pub_counts()
    lambert   = _load_lambert_tfs()
    perturbed = _load_perturbed_tfs()
    sub_gene_data = []
    for g in genes:
        if g in perturbed:
            tf_status = 'Active TF'
        elif g in lambert:
            tf_status = 'TF'
        else:
            tf_status = 'Gene'
        sub_gene_data.append({'name': g, 'pubs': pub_counts.get(g, 0), 'tf_status': tf_status})
    raw_neighbors = rows_to_dicts(db.execute(
        "SELECT source, target, mean_z, n_pairs FROM submodule_edges "
        "WHERE source=? OR target=? ORDER BY ABS(mean_z) DESC LIMIT 40", (name, name)))
    for e in raw_neighbors:
        e["other"] = e["target"] if e["source"] == name else e["source"]
    expr = rows_to_dicts(db.execute(
        f"SELECT timepoint, mean_tpm, sd_tpm, n_genes FROM submodule_expression WHERE submodule=? {_tp_order()}",
        (name,)))
    timepoints = [r['timepoint'] for r in expr]
    expr_z = _compute_mean_zscore_profile(db, genes, timepoints)
    highlight_terms = [r[0] for r in db.execute(
        "SELECT term_name FROM module_enrichment "
        "WHERE module=? AND module_type='submodule' AND significant='TRUE' "
        "AND source IN ('GO:BP','GO:CC','GO:MF') "
        "ORDER BY p_value LIMIT 4", (name,)).fetchall()]
    tf_rows = db.execute(
        "SELECT tf, direction FROM perturbation_edges "
        "WHERE module=? AND module_type='submodule' "
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
    # Use LLM-generated description if available, otherwise fall back to auto-generated
    llm_descs = _load_submodule_descriptions()
    llm_desc = llm_descs.get(name)
    top_neighbor_ids = [e["other"] for e in raw_neighbors[:10]]
    if top_neighbor_ids:
        placeholders = ','.join('?' * len(top_neighbor_ids))
        supermod_rows = {r['id']: r['supermodule'] for r in db.execute(
            f"SELECT id, supermodule FROM submodule_nodes WHERE id IN ({placeholders})",
            top_neighbor_ids).fetchall()}
    else:
        supermod_rows = {}
    similar_submodules = [
        {'name': sid, 'color': MODULE_COLORS.get(supermod_rows.get(sid, ''), '#4a56b0')}
        for sid in top_neighbor_ids
    ]
    desc = {
        'title': highlight_terms[0] if highlight_terms else name,
        'summary': llm_desc or (' '.join(summary_parts) or None),
        'highlight_terms': highlight_terms,
        'highlight_genes': notable,
        'top_tfs': top_tfs_auto,
        'similar_submodules': similar_submodules,
    }
    tfs, modules = _nav_lists(db)
    supermodule = node.get('supermodule', '')
    card_color = MODULE_COLORS.get(supermodule, node.get('color', '#4a56b0'))
    return render_template("perturbseq/module.html",
        module_type='sub', name=name, color=node.get('color', '#4a56b0'), card_color=card_color,
        desc=desc, expr=expr, expr_z=expr_z, profile=None, parent=supermodule, node=node,
        submodules=[], neighbors=raw_neighbors, genes=genes,
        sub_gene_data=sub_gene_data,
        notable_genes=notable,
        top_tfs=top_tfs_auto, assoc=[], assoc_label='', enrichment=[],
        tfs=tfs, modules=modules)


def _gc_page(name):
    internal_id = 'gene_cluster_' + name[2:]
    if internal_id in HIDDEN_TC:
        return render_template("perturbseq/404.html", message=f"Cluster not found: {name}"), 404
    tc = get_tc_db()
    profile = rows_to_dicts(tc.execute(
        "SELECT timepoint, value FROM gene_cluster_profiles WHERE cluster=? ORDER BY rowid",
        (internal_id,)).fetchall())
    if not profile:
        return render_template("perturbseq/404.html", message=f"Cluster not found: {name}"), 404
    pert_tfs = {r['tf']: r for r in rows_to_dicts(tc.execute(
        "SELECT tf, nes, padj FROM perturbation_edges WHERE gene_cluster=?",
        (internal_id,)).fetchall())}
    bind_tfs = {
        r[0]: round(r[1], 3)
        for r in tc.execute("""
            SELECT b.tf, MAX(b.odds_ratio) as odds_ratio
            FROM binding_edges b
            JOIN peak_gene_edges p ON b.peak_cluster = p.peak_cluster
            WHERE p.gene_cluster = ? AND b.fdr < 0.05 AND p.fdr < 0.05
            GROUP BY b.tf
        """, (internal_id,)).fetchall()
    }
    top_tfs = []
    for tf, p in pert_tfs.items():
        b = bind_tfs.get(tf)
        top_tfs.append({
            'tf': tf, 'nes': p['nes'], 'padj': p['padj'],
            'binding': '✓' if b is not None else '',
            'odds_ratio': b,
            'evidence': 'both' if b is not None else 'perturbation',
        })
    for tf, or_val in bind_tfs.items():
        if tf not in pert_tfs:
            top_tfs.append({
                'tf': tf, 'nes': None, 'padj': None,
                'binding': '✓', 'odds_ratio': or_val, 'evidence': 'binding',
            })
    top_tfs.sort(key=lambda r: abs(r['nes'] or 0), reverse=True)
    assoc_raw = rows_to_dicts(tc.execute(
        "SELECT peak_cluster, odds_ratio, fdr FROM peak_gene_edges WHERE gene_cluster=? ORDER BY odds_ratio DESC",
        (internal_id,)).fetchall())
    assoc = [{'name': _tc_display(r['peak_cluster']),
              'url': url_for('perturbseq.module_page', name=_tc_display(r['peak_cluster'])),
              'odds_ratio': r['odds_ratio'], 'fdr': r['fdr']}
             for r in assoc_raw if r['peak_cluster'] not in HIDDEN_TC]
    rev_map = {v: k for k, v in _load_gene_name_map().items()}
    members_ensg = [r[0] for r in tc.execute(
        "SELECT gene FROM gene_cluster_members WHERE cluster=? ORDER BY membership DESC LIMIT 500",
        (internal_id,)).fetchall()]
    genes = [rev_map.get(e, e) for e in members_ensg]
    db = get_db()
    desc_row = db.execute("SELECT * FROM module_descriptions WHERE module=?", (name,)).fetchone()
    desc = None
    if desc_row:
        desc = dict(desc_row)
        try:
            desc["highlight_terms"] = json.loads(desc["highlight_terms"])
        except (json.JSONDecodeError, TypeError):
            desc["highlight_terms"] = []
        try:
            desc["top_tfs"] = json.loads(desc["top_tfs"])
        except (json.JSONDecodeError, TypeError):
            desc["top_tfs"] = []
        desc["highlight_genes"] = _top_genes_by_pub(genes, n=5)
    pub_counts = _load_gene_pub_counts()
    lambert    = _load_lambert_tfs()
    perturbed  = _load_perturbed_tfs()
    gc_gene_data = []
    for g in genes:
        if g.startswith('ENSG'):
            continue
        if g in perturbed:
            tf_status = 'Active TF'
        elif g in lambert:
            tf_status = 'TF'
        else:
            tf_status = 'Gene'
        gc_gene_data.append({'name': g, 'pubs': pub_counts.get(g, 0), 'tf_status': tf_status})
    timepoints = [r['timepoint'] for r in profile]
    named_genes = [g for g in genes if not g.startswith('ENSG')]
    expr_z = _compute_mean_zscore_profile(db, named_genes, timepoints)
    expr = _compute_mean_tpm_profile(db, named_genes, timepoints)
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
    tc = get_tc_db()
    profile = rows_to_dicts(tc.execute(
        "SELECT timepoint, value FROM peak_cluster_profiles WHERE cluster=? ORDER BY rowid",
        (internal_id,)).fetchall())
    if not profile:
        return render_template("perturbseq/404.html", message=f"Cluster not found: {name}"), 404
    top_chipseq = rows_to_dicts(tc.execute(
        "SELECT tf, odds_ratio, fdr, source FROM binding_edges "
        "WHERE peak_cluster=? AND source='chipseq' ORDER BY odds_ratio DESC LIMIT 50",
        (internal_id,)).fetchall())
    top_footprint = rows_to_dicts(tc.execute(
        "SELECT tf, odds_ratio, fdr, source FROM binding_edges "
        "WHERE peak_cluster=? AND source='footprint' ORDER BY odds_ratio DESC LIMIT 50",
        (internal_id,)).fetchall())
    assoc_raw = rows_to_dicts(tc.execute(
        "SELECT gene_cluster, odds_ratio, fdr FROM peak_gene_edges WHERE peak_cluster=? ORDER BY odds_ratio DESC",
        (internal_id,)).fetchall())
    assoc = [{'name': _tc_display(r['gene_cluster']),
              'url': url_for('perturbseq.module_page', name=_tc_display(r['gene_cluster'])),
              'odds_ratio': r['odds_ratio'], 'fdr': r['fdr']}
             for r in assoc_raw if r['gene_cluster'] not in HIDDEN_TC]
    members = [r[0] for r in tc.execute(
        "SELECT peak FROM peak_cluster_members WHERE cluster=? ORDER BY membership DESC LIMIT 100",
        (internal_id,)).fetchall()]
    db = get_db()
    tfs, modules = _nav_lists(db)
    return render_template("perturbseq/module.html",
        module_type='pc', name=name, color='#e07b39', card_color='#e07b39',
        desc=None, expr=None, expr_z=[], profile=profile, parent=None, node=None,
        submodules=[], neighbors=[], genes=members,
        notable_genes=[],
        top_tfs=top_chipseq + top_footprint, assoc=assoc, assoc_label='Gene Clusters',
        enrichment=[], tfs=tfs, modules=modules)


@perturbseq_bp.route("/go/<path:term_name>")
def go_page(term_name):
    return _go_page(term_name)


def _go_page(term_name):
    import re as _re
    db = get_db()

    # If given a bare GO ID (e.g. "GO:0007015"), resolve to full term name
    if _re.match(r'^GO:\d+$', term_name):
        row = db.execute(
            "SELECT DISTINCT term FROM go_genes WHERE term LIKE ?",
            ('%(' + term_name + ')%',)).fetchone()
        if row:
            return redirect(url_for('perturbseq.go_page', term_name=row[0]), 301)
        return render_template("perturbseq/404.html",
                               message=f"GO term not found: {term_name}"), 404

    # Get member genes for this GO term
    genes = [r[0] for r in db.execute(
        "SELECT gene FROM go_genes WHERE term=? ORDER BY gene", (term_name,))]
    if not genes:
        return render_template("perturbseq/404.html",
                               message=f"GO term not found: {term_name}"), 404

    # Extract GO ID from term name, e.g. "Actin Filament Organization (GO:0007015)" → "GO:0007015"
    go_id_match = _re.search(r'\((GO:\d+)\)', term_name)
    go_id = go_id_match.group(1) if go_id_match else None

    # Find modules enriched for this GO term (reverse lookup via term_id)
    enriched_modules = []
    if go_id:
        enriched_modules = rows_to_dicts(db.execute(
            "SELECT module, module_type, p_value, term_size, intersection_size, "
            "precision_val, recall, query_size FROM module_enrichment "
            "WHERE term_id=? AND source='GO:BP' ORDER BY p_value",
            (go_id,)))
        for row in enriched_modules:
            mt = row['module_type']
            if mt == 'supermodule':
                row['color'] = MODULE_COLORS.get(row['module'], '#4a56b0')
                row['url'] = url_for('perturbseq.module_page', name=row['module'])
            elif mt == 'submodule':
                sup = row['module'].rsplit('.', 1)[0] if '.' in row['module'] else ''
                row['color'] = MODULE_COLORS.get(sup, '#4a56b0')
                row['url'] = url_for('perturbseq.module_page', name=row['module'])
            elif mt == 'timecourse':
                row['color'] = '#4e79a7'
                row['url'] = url_for('perturbseq.module_page', name=_tc_display(row['module']))
                row['module'] = _tc_display(row['module'])
            else:
                row['color'] = '#4a56b0'
                row['url'] = '#'

    # TFs that perturb this GO term gene set
    pert_rows = rows_to_dicts(db.execute(
        "SELECT tf, mean_NES, n_sig_gRNA, padj FROM perturbation_edges "
        "WHERE module=? AND module_type='go' ORDER BY ABS(mean_NES) DESC",
        (term_name,)))

    # TFs that bind this GO term gene set
    bind_map = {r['tf']: r for r in rows_to_dicts(db.execute(
        "SELECT tf, odds_ratio FROM binding_edges "
        "WHERE module=? AND module_type='go'", (term_name,)))}

    # Merge perturbation + binding into a unified TF table
    tf_rows = _merge_pert_bind_edges(pert_rows, bind_map)

    # Gene table data: pub counts + TF status
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
    """Return ATAC peaks bound by tf_name that also link to at least one gene region."""
    rows = db.execute("""
        SELECT DISTINCT
            ap.atac_peak_id, ap.chr, ap.start, ap.end,
            gr.gene_name
        FROM tf_datasets    td
        JOIN tf_peaks          tp  ON td.dataset_id     = tp.dataset_id
        JOIN atac_tf_overlaps  ato ON tp.peak_id         = ato.peak_id
        JOIN atac_peaks        ap  ON ato.atac_peak_id   = ap.atac_peak_id
        JOIN tf_gene_overlaps  tgo ON tp.peak_id         = tgo.peak_id
        JOIN gene_regions      gr  ON tgo.gene_region_id = gr.gene_region_id
        WHERE td.tf_gene_name = ?
        ORDER BY ap.chr, ap.start, gr.gene_name
    """, (tf_name,)).fetchall()
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
    import json as _json

    # Expression profile
    expr = rows_to_dicts(db.execute(
        f"SELECT timepoint, mean_tpm, replicates FROM gene_expression WHERE gene=? {_tp_order()}",
        (gene_name,)))
    if not expr:
        return render_template("perturbseq/404.html",
                               message=f"Gene not found: {gene_name}"), 404
    for row in expr:
        try:
            row['replicates'] = _json.loads(row['replicates'])
        except (json.JSONDecodeError, TypeError):
            row['replicates'] = []

    # Module membership (as a gene)
    _mem = db.execute(
        "SELECT submodule, supermodule FROM submodule_genes WHERE gene=?",
        (gene_name,)).fetchone()
    if _mem and _mem['supermodule'] == 'unassigned':
        _mod = db.execute(
            "SELECT module FROM module_genes WHERE gene=?", (gene_name,)).fetchone()
        membership = {'submodule': None, 'supermodule': _mod['module']} if _mod else None
    else:
        membership = dict(_mem) if _mem else None

    # Developmental gene cluster membership (timecourse)
    tc_clusters = []
    gene_map = _load_gene_name_map()
    ensg_id = gene_map.get(gene_name)
    if ensg_id:
        try:
            tc_db = get_tc_db()
            for r in tc_db.execute(
                "SELECT cluster, membership FROM gene_cluster_members WHERE gene=? "
                "ORDER BY membership DESC", (ensg_id,)).fetchall():
                tc_clusters.append({
                    "cluster": r[0],
                    "label": r[0].replace("gene_cluster_", "GC"),
                    "membership": r[1],
                    "color": TC_CLUSTER_COLORS.get(r[0], "#8090b0"),
                    "url": url_for('perturbseq.module_page', name=_tc_display(r[0])),
                })
        except Exception:
            pass

    # TF-specific data (only if this gene is an active TF)
    is_tf = bool(db.execute(
        "SELECT 1 FROM perturbation_edges WHERE tf=? "
        "UNION ALL SELECT 1 FROM binding_edges WHERE tf=? LIMIT 1",
        (gene_name, gene_name)).fetchone())
    tf_desc_row = db.execute(
        "SELECT name, summary FROM tf_descriptions WHERE tf=?", (gene_name,)).fetchone()
    tf_desc = dict(tf_desc_row) if tf_desc_row else None

    tf_reg_modules = []
    submodules_by_mod = {}
    pert_by_mod = {}
    tc_reg_clusters = []

    if is_tf:
        pert_rows = rows_to_dicts(db.execute(
            "SELECT module, mean_NES, n_sig_gRNA, padj FROM perturbation_edges "
            "WHERE tf=? AND module_type='supermodule'",
            (gene_name,)))
        bind_by_mod = {r['module']: r for r in rows_to_dicts(db.execute(
            "SELECT module, odds_ratio FROM binding_edges "
            "WHERE tf=? AND module_type='supermodule'", (gene_name,)))}
        pert_mods = {r['module'] for r in pert_rows}
        pert_by_mod = {r['module']: r['mean_NES'] for r in pert_rows}

        tf_reg_modules = _merge_pert_bind_edges(pert_rows, bind_by_mod, id_key='module')
        for row in tf_reg_modules:
            row['color'] = MODULE_COLORS.get(row['module'], '#4a56b0')
        tf_reg_modules.sort(key=lambda r: abs(r['mean_NES'] or 0), reverse=True)

        for mod in list(pert_mods | set(bind_by_mod.keys())):
            subs = rows_to_dicts(db.execute(
                "SELECT id, n_genes, color, within_mean_z FROM submodule_nodes "
                "WHERE supermodule=? ORDER BY id", (mod,)))
            if subs:
                sup_color = MODULE_COLORS.get(mod, '#4a56b0')
                for s in subs:
                    s['supermodule_color'] = sup_color
                submodules_by_mod[mod] = subs

        try:
            tc_db = get_tc_db()
            HIDDEN = {'gene_cluster_7'}
            pert_gc = {
                r[0]: {'nes': r[1], 'padj': r[2]}
                for r in tc_db.execute(
                    "SELECT gene_cluster, nes, padj FROM perturbation_edges WHERE tf=?",
                    (gene_name,)).fetchall()
                if r[0] not in HIDDEN
            }
            bind_gc = {
                r[0]: round(r[1], 3)
                for r in tc_db.execute("""
                    SELECT p.gene_cluster, MAX(b.odds_ratio) as odds_ratio
                    FROM binding_edges b
                    JOIN peak_gene_edges p ON b.peak_cluster = p.peak_cluster
                    WHERE b.tf = ? AND b.fdr < 0.05 AND p.fdr < 0.05
                    GROUP BY p.gene_cluster
                """, (gene_name,)).fetchall()
                if r[0] not in HIDDEN
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
                    "odds_ratio": b,
                    "evidence": evidence,
                })
            tc_reg_clusters.sort(key=lambda r: abs(r['nes'] or 0), reverse=True)
        except Exception:
            pass

    card_color = MODULE_COLORS.get(
        membership['supermodule'] if membership else '', '#4a56b0')
    submodules_flat = []
    if is_tf:
        # Gather submodule-level perturbation + binding for this TF
        sub_pert = {r['module']: r for r in rows_to_dicts(db.execute(
            "SELECT module, mean_NES, n_sig_gRNA, padj FROM perturbation_edges "
            "WHERE tf=? AND module_type='submodule'", (gene_name,)))}
        sub_bind = {r['module']: r for r in rows_to_dicts(db.execute(
            "SELECT module, odds_ratio FROM binding_edges "
            "WHERE tf=? AND module_type='submodule'", (gene_name,)))}
        all_sub_ids = set(sub_pert.keys()) | set(sub_bind.keys())
        # Get submodule metadata (n_genes, supermodule)
        if all_sub_ids:
            placeholders = ','.join('?' * len(all_sub_ids))
            sub_meta = {r['id']: r for r in rows_to_dicts(db.execute(
                f"SELECT id, supermodule, n_genes, color FROM submodule_nodes "
                f"WHERE id IN ({placeholders})", list(all_sub_ids)))}
        else:
            sub_meta = {}
        for sub_id in sorted(all_sub_ids):
            meta = sub_meta.get(sub_id, {})
            sup = meta.get('supermodule', sub_id.rsplit('.', 1)[0] if '.' in sub_id else '')
            sup_color = MODULE_COLORS.get(sup, '#4a56b0')
            p = sub_pert.get(sub_id)
            b = sub_bind.get(sub_id)
            if p and b:
                evidence = 'both'
            elif p:
                evidence = 'perturbation'
            else:
                evidence = 'binding'
            direction = None
            if p:
                direction = 'Up' if p['mean_NES'] > 0 else 'Down'
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
                'evidence': evidence,
            })
        submodules_flat.sort(key=lambda x: abs(x['mean_NES'] or 0), reverse=True)
    tfs     = [r[0] for r in db.execute("SELECT DISTINCT tf FROM tf_descriptions ORDER BY 1")]
    modules = [r[0] for r in db.execute("SELECT DISTINCT module FROM module_descriptions ORDER BY 1")]
    reg_elements      = _query_gene_elements(db, gene_name)
    tf_element_count  = _count_tf_elements(db, gene_name) if is_tf else 0
    _locus_row = db.execute(
        "SELECT chr, MIN(start) AS locus_start, MAX(end) AS locus_end "
        "FROM gene_regions WHERE gene_name=? AND gene_region_subtype IN ('proximal','distal') "
        "GROUP BY chr ORDER BY (MAX(end)-MIN(start)) DESC LIMIT 1",
        (gene_name,)).fetchone()
    gene_locus = dict(_locus_row) if _locus_row else None
    if gene_locus:
        gene_locus['chr'] = _norm_chr(gene_locus['chr'])
    return render_template("perturbseq/gene.html",
        gene=gene_name, expr=expr,
        membership=membership,
        tc_clusters=tc_clusters,
        is_tf=is_tf, tf_desc=tf_desc,
        card_color=card_color,
        tf_reg_modules=tf_reg_modules,
        submodules_by_mod=submodules_by_mod, pert_by_mod=pert_by_mod,
        submodules_flat=submodules_flat,
        tc_reg_clusters=tc_reg_clusters,
        tfs=tfs, modules=modules,
        reg_elements=reg_elements,
        tf_element_count=tf_element_count,
        gene_locus=gene_locus)


def _count_tf_elements(db, tf_name: str) -> int:
    # CTCF has 26 M rows in tf_peaks; any full-scan query times out.
    # Use LIMIT 502 so SQLite's hash-DISTINCT terminates early after finding 502
    # distinct atac_peak_ids.  Returns exact count ≤501, or -1 as a sentinel
    # meaning "more than 501 — use server-side table mode without total count."
    rows = db.execute("""
        SELECT DISTINCT ato.atac_peak_id
        FROM tf_datasets      td
        JOIN tf_peaks         tp  ON td.dataset_id = tp.dataset_id
        JOIN atac_tf_overlaps ato ON tp.peak_id    = ato.peak_id
        WHERE td.tf_gene_name = ?
          AND EXISTS (SELECT 1 FROM tf_gene_overlaps WHERE peak_id = tp.peak_id)
        LIMIT 502
    """, (tf_name,)).fetchall()
    return len(rows) if len(rows) < 502 else -1


def _query_tf_elements_paged(db, tf_name: str, start: int, length: int,
                              q: str, order_col: int, order_dir: str) -> dict:
    """Paginated + filtered TF elements."""
    safe_dir = 'DESC' if order_dir.upper() == 'DESC' else 'ASC'
    order_map = {
        0: f'ap.chr {safe_dir}, ap.start {safe_dir}',
        1: f'(ap.end - ap.start) {safe_dir}',
    }
    order_by = order_map.get(order_col, f'ap.chr {safe_dir}, ap.start {safe_dir}')

    if q:
        # Search path: need to join gene_regions to filter by gene name.
        # Kept as the original join but only used when the user is actively searching.
        like = f'%{q}%'
        base_from = """
            FROM tf_datasets      td
            JOIN tf_peaks         tp  ON td.dataset_id     = tp.dataset_id
            JOIN atac_tf_overlaps ato ON tp.peak_id        = ato.peak_id
            JOIN atac_peaks       ap  ON ato.atac_peak_id  = ap.atac_peak_id
            JOIN tf_gene_overlaps tgo ON tp.peak_id        = tgo.peak_id
            JOIN gene_regions     gr  ON tgo.gene_region_id = gr.gene_region_id
            WHERE td.tf_gene_name = ?
            GROUP BY ap.atac_peak_id
            HAVING (ap.chr LIKE ? OR GROUP_CONCAT(DISTINCT gr.gene_name) LIKE ?)
        """
        params = [tf_name, like, like]
        count_row = db.execute(
            f"SELECT COUNT(*) FROM (SELECT ap.atac_peak_id {base_from})", params
        ).fetchone()
        filtered_count = count_row[0] if count_row else 0
        rows = db.execute(f"""
            SELECT ap.atac_peak_id, ap.chr, ap.start, ap.end,
                   GROUP_CONCAT(DISTINCT gr.gene_name) AS genes_str
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

    # Fast path (no search): LIMIT-502 capped count + no ORDER BY so SQLite's
    # hash-DISTINCT can early-terminate.  ORDER BY forces a full sort of all
    # 28k+ elements for genome-wide binders like CTCF (>30 s timeout).
    count_rows = db.execute("""
        SELECT DISTINCT ato.atac_peak_id
        FROM tf_datasets      td
        JOIN tf_peaks         tp  ON td.dataset_id = tp.dataset_id
        JOIN atac_tf_overlaps ato ON tp.peak_id    = ato.peak_id
        WHERE td.tf_gene_name = ?
          AND EXISTS (SELECT 1 FROM tf_gene_overlaps WHERE peak_id = tp.peak_id)
        LIMIT 502
    """, [tf_name]).fetchall()
    filtered_count = len(count_rows) if len(count_rows) < 502 else -1

    # ORDER BY only if the user explicitly requested a sort column; otherwise
    # drop it so LIMIT/OFFSET can early-terminate on the hash-DISTINCT.
    explicit_order = order_map.get(order_col) if order_col != 0 else None
    order_clause   = f"ORDER BY {explicit_order}" if explicit_order else ""

    page_rows = db.execute(f"""
        SELECT DISTINCT ap.atac_peak_id, ap.chr, ap.start, ap.end
        FROM tf_datasets      td
        JOIN tf_peaks         tp  ON td.dataset_id    = tp.dataset_id
        JOIN atac_tf_overlaps ato ON tp.peak_id       = ato.peak_id
        JOIN atac_peaks       ap  ON ato.atac_peak_id = ap.atac_peak_id
        WHERE td.tf_gene_name = ?
          AND EXISTS (SELECT 1 FROM tf_gene_overlaps WHERE peak_id = tp.peak_id)
        {order_clause}
        LIMIT ? OFFSET ?
    """, [tf_name, length, start]).fetchall()

    if not page_rows:
        return {'data': [], 'filtered': filtered_count}

    # Fetch linked gene names for only the page's elements.
    # Start from atac_tf_overlaps filtered to the page's peak_ids — this uses
    # idx_atac_tf_overlaps_atac_peak_id and avoids scanning all 26M CTCF tf_peaks.
    peak_ids = [r['atac_peak_id'] for r in page_rows]
    placeholders = ','.join('?' * len(peak_ids))
    gene_rows = db.execute(f"""
        SELECT DISTINCT ato.atac_peak_id, gr.gene_name
        FROM atac_tf_overlaps ato
        JOIN tf_peaks         tp  ON ato.peak_id        = tp.peak_id
        JOIN tf_datasets      td  ON tp.dataset_id      = td.dataset_id
        JOIN tf_gene_overlaps tgo ON ato.peak_id        = tgo.peak_id
        JOIN gene_regions     gr  ON tgo.gene_region_id = gr.gene_region_id
        WHERE ato.atac_peak_id IN ({placeholders})
          AND td.tf_gene_name = ?
    """, peak_ids + [tf_name]).fetchall()
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


@perturbseq_bp.route("/api/gene/<gene_name>")
def api_gene(gene_name):
    db = get_db()
    import json as _json
    expr = rows_to_dicts(db.execute(
        f"SELECT timepoint, mean_tpm, replicates FROM gene_expression WHERE gene=? {_tp_order()}",
        (gene_name,)))
    if not expr:
        return jsonify({"error": "not found"}), 404
    for row in expr:
        try:
            row['replicates'] = _json.loads(row['replicates'])
        except (json.JSONDecodeError, TypeError):
            row['replicates'] = []
    membership = db.execute(
        "SELECT submodule, supermodule FROM submodule_genes WHERE gene=?",
        (gene_name,)).fetchone()
    return jsonify({"gene": gene_name, "expr": expr,
                    "membership": dict(membership) if membership else None})


@perturbseq_bp.route("/api/submodule/<name>/gene_network")
def api_submodule_gene_network(name):
    db = get_db()
    edges = rows_to_dicts(db.execute(
        "SELECT gene1, gene2, z_score FROM submodule_gene_edges "
        "WHERE submodule=? ORDER BY z_score DESC", (name,)))
    if not edges:
        return jsonify({"edges": [], "z_min": 0, "z_max": 0, "z_p50": 0})
    zscores = [e["z_score"] for e in edges]
    zscores_sorted = sorted(zscores)
    n = len(zscores_sorted)
    return jsonify({
        "edges": edges,
        "z_min": zscores_sorted[0],
        "z_max": zscores_sorted[-1],
        "z_p50": zscores_sorted[n // 2],
        "z_p75": zscores_sorted[int(n * 0.75)],
    })


@perturbseq_bp.route("/api/module/<mod_name>/submodule_network")
def api_module_submodule_network(mod_name):
    db = get_db()
    nodes_raw = rows_to_dicts(db.execute(
        "SELECT id, n_genes, color, within_mean_z FROM submodule_nodes WHERE supermodule=?",
        (mod_name,)))
    if not nodes_raw:
        return jsonify({"nodes": [], "edges": []})
    node_ids = {n["id"] for n in nodes_raw}
    all_edges = rows_to_dicts(db.execute(
        "SELECT source, target, mean_z, n_pairs FROM submodule_edges "
        "WHERE source LIKE ? OR target LIKE ?",
        (mod_name + '.%', mod_name + '.%')))
    edges_filtered = [e for e in all_edges
                      if e["source"] in node_ids and e["target"] in node_ids]
    return jsonify({"nodes": nodes_raw, "edges": edges_filtered})


@perturbseq_bp.route("/api/submodule/<name>")
def api_submodule(name):
    db = get_db()
    node = db.execute("SELECT * FROM submodule_nodes WHERE id=?", (name,)).fetchone()
    if not node:
        return jsonify({"error": "not found"}), 404
    genes = [r[0] for r in db.execute(
        "SELECT gene FROM submodule_genes WHERE submodule=? ORDER BY gene", (name,))]
    neighbors = rows_to_dicts(db.execute(
        "SELECT source, target, mean_z, n_pairs FROM submodule_edges "
        "WHERE source=? OR target=? ORDER BY ABS(mean_z) DESC LIMIT 40",
        (name, name)))
    for e in neighbors:
        e["other"] = e["target"] if e["source"] == name else e["source"]
    return jsonify({"node": dict(node), "genes": genes, "neighbors": neighbors})


@perturbseq_bp.route("/api/edges")
def api_edges():
    db = get_db()
    min_strength = float(request.args.get("min_strength", 0))
    evidence = request.args.get("evidence", "all")
    tf_filter = request.args.get("tf", "").upper()
    mod_filter = request.args.get("module", "").upper()

    expressed = _load_expressed_tfs()

    pert = rows_to_dicts(db.execute(
        "SELECT tf, module, mean_NES, n_sig_gRNA, padj FROM perturbation_edges").fetchall())
    bind_map = {}
    for r in db.execute("SELECT tf, module, odds_ratio, padj FROM binding_edges").fetchall():
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
        "SELECT module, mean_NES, n_sig_gRNA, padj FROM perturbation_edges WHERE tf=?",
        (tf_name,)).fetchall())
    bind = rows_to_dicts(db.execute(
        "SELECT module, odds_ratio, log2or, padj FROM binding_edges WHERE tf=?",
        (tf_name,)).fetchall())
    return jsonify({"perturbation": pert, "binding": bind})


@perturbseq_bp.route("/api/module/<mod_name>")
def api_module(mod_name):
    db = get_db()
    pert = rows_to_dicts(db.execute(
        "SELECT tf, mean_NES, n_sig_gRNA, padj FROM perturbation_edges WHERE module=?",
        (mod_name,)).fetchall())
    bind = rows_to_dicts(db.execute(
        "SELECT tf, odds_ratio, log2or, padj FROM binding_edges WHERE module=?",
        (mod_name,)).fetchall())
    return jsonify({"perturbation": pert, "binding": bind})


@perturbseq_bp.route("/api/module/<mod_name>/detail")
def api_module_detail(mod_name):
    db = get_db()
    enrichment = rows_to_dicts(db.execute(
        "SELECT term_id, source, term_name, p_value, term_size, intersection_size, precision_val, recall "
        "FROM module_enrichment WHERE module=? AND source IN ('GO:BP','GO:CC','GO:MF')", (mod_name,)).fetchall())
    genes = [r[0] for r in db.execute(
        "SELECT gene FROM module_genes WHERE module=? ORDER BY gene", (mod_name,)).fetchall()]
    pub_counts = _load_gene_pub_counts()
    gene_pubs = {g: pub_counts.get(g, 0) for g in genes}
    lambert   = _load_lambert_tfs()
    perturbed = _load_perturbed_tfs()
    gene_tf_status = {}
    for g in genes:
        if g in perturbed:
            gene_tf_status[g] = 'Active TF'
        elif g in lambert:
            gene_tf_status[g] = 'TF'
        else:
            gene_tf_status[g] = 'Gene'
    sub_colors = {r[0]: r[1] for r in db.execute(
        "SELECT id, color FROM submodule_nodes WHERE supermodule=?", (mod_name,)).fetchall()}
    gene_submodule = {r[0]: {"id": r[1], "color": sub_colors.get(r[1], "#888")}
                      for r in db.execute(
        "SELECT gene, submodule FROM submodule_genes WHERE submodule IN "
        "(SELECT id FROM submodule_nodes WHERE supermodule=?)", (mod_name,)).fetchall()}
    return jsonify({"enrichment": enrichment, "genes": genes,
                    "pub_counts": gene_pubs, "gene_tf_status": gene_tf_status,
                    "gene_submodule": gene_submodule})


@perturbseq_bp.route("/api/module/<mod_name>/genes")
def api_module_genes(mod_name):
    db = get_db()
    genes = [r[0] for r in db.execute(
        "SELECT gene FROM module_genes WHERE module=? ORDER BY gene", (mod_name,)).fetchall()]
    return jsonify(genes)


@perturbseq_bp.route("/api/module/<mod_name>/enrichment")
def api_module_enrichment(mod_name):
    db = get_db()
    rows = rows_to_dicts(db.execute(
        "SELECT term_id, term_name, source, p_value, term_size, intersection_size, precision_val, recall "
        "FROM module_enrichment WHERE module=? AND source IN ('GO:BP','GO:CC','GO:MF') ORDER BY p_value", (mod_name,)).fetchall())
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
                "SELECT DISTINCT module FROM module_descriptions WHERE module LIKE ? ORDER BY module LIMIT 20",
                (q + "%",)).fetchall()
            return jsonify([r[0] for r in rows])

        results = []

        # 1. Active TFs — in perturbseq dataset, have a dedicated TF page
        tf_rows = db.execute(
            "SELECT DISTINCT tf FROM tf_descriptions WHERE tf LIKE ? ORDER BY tf LIMIT 8",
            (q + "%",)).fetchall()
        active_tf_set = {r[0] for r in tf_rows}
        results.extend({"name": r[0], "type": "active_tf"} for r in tf_rows)

        # 2. Modules
        mod_rows = db.execute(
            "SELECT DISTINCT module FROM module_descriptions WHERE module LIKE ? ORDER BY module LIMIT 5",
            (q + "%",)).fetchall()
        results.extend({"name": r[0], "type": "module"} for r in mod_rows)

        # 2b. Submodules
        sub_rows = db.execute(
            "SELECT DISTINCT id FROM submodule_nodes WHERE id LIKE ? ORDER BY id LIMIT 5",
            (q + "%",)).fetchall()
        results.extend({"name": r[0], "type": "submodule"} for r in sub_rows)

        # 2d. GO terms
        go_rows = db.execute(
            "SELECT DISTINCT term FROM go_genes WHERE term LIKE ? ORDER BY LENGTH(term) LIMIT 5",
            ('%' + q + '%',)).fetchall()
        results.extend({"name": r[0], "type": "go_term",
                        "url": "/perturbseq/go/" + r[0]} for r in go_rows)

        # 2c. TC gene / peak clusters
        try:
            tc = get_tc_db()
            q_up = q.upper()
            for cid in [r[0] for r in tc.execute(
                    "SELECT DISTINCT cluster FROM gene_cluster_profiles ORDER BY cluster")]:
                if cid in HIDDEN_TC: continue
                dn = _tc_display(cid)
                if dn.upper().startswith(q_up):
                    results.append({"name": dn, "type": "gc_cluster",
                                    "url": "/perturbseq/module/" + dn})
            for cid in [r[0] for r in tc.execute(
                    "SELECT DISTINCT cluster FROM peak_cluster_profiles ORDER BY cluster")]:
                if cid in HIDDEN_TC: continue
                dn = _tc_display(cid)
                if dn.upper().startswith(q_up):
                    results.append({"name": dn, "type": "pc_cluster",
                                    "url": "/perturbseq/module/" + dn})
        except Exception:
            pass

        # 3. Lambert TFs — known TFs not in the perturbseq dataset
        all_lambert = _load_lambert_tfs()
        q_upper = q.upper()
        lambert_matches = sorted(
            name for name in all_lambert
            if name.upper().startswith(q_upper) and name not in active_tf_set
        )[:8]
        results.extend({"name": name, "type": "tf"} for name in lambert_matches)

        # 4. Non-TF genes
        gene_rows = db.execute(
            "SELECT DISTINCT gene FROM gene_expression WHERE gene LIKE ? ORDER BY gene LIMIT 20",
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

def _subtype_to_link_type(subtype: str, hicar_support: str) -> str:
    if subtype == 'proximal':
        return 'tss_proximal_1kb'
    if subtype == 'distal':
        return 'tss_distal_10kb'
    if hicar_support and hicar_support != 'no_HiCAR_support':
        return 'distal_multiome_hicar'
    return 'distal_multiome'


def _norm_chr(chrom: str) -> str:
    """Ensure chromosome has 'chr' prefix (raw data uses bare numbers)."""
    if chrom and not chrom.startswith('chr'):
        return 'chr' + chrom
    return chrom


def _query_tf_gene_link(db, tf: str, gene: str) -> list:
    """Return list of regulatory element dicts for a (TF, gene) pair."""
    rows = db.execute("""
        SELECT
            ap.atac_peak_id, ap.chr, ap.start, ap.end,
            tp.source  AS peak_source,
            td.dataset_id, td.dataset, td.cell_type,
            gr.gene_region_subtype, gr.multiome_hicar_support,
            gr.start AS gr_start, gr.end AS gr_end
        FROM tf_datasets td
        JOIN tf_peaks          tp  ON td.dataset_id      = tp.dataset_id
        JOIN tf_gene_overlaps  tgo ON tp.peak_id         = tgo.peak_id
        JOIN gene_regions      gr  ON tgo.gene_region_id = gr.gene_region_id
        JOIN atac_tf_overlaps  ato ON tp.peak_id         = ato.peak_id
        JOIN atac_peaks        ap  ON ato.atac_peak_id   = ap.atac_peak_id
        WHERE td.tf_gene_name = ? AND gr.gene_name = ?
    """, (tf, gene)).fetchall()

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

        link_key = (r['gene_region_subtype'], r['multiome_hicar_support'])
        if link_key not in el['_seen_links']:
            el['_seen_links'].add(link_key)
            link_type = _subtype_to_link_type(r['gene_region_subtype'], r['multiome_hicar_support'])
            # distance_to_tss: midpoint of ATAC peak relative to midpoint of gene region
            atac_mid = (r['start'] + r['end']) / 2
            gr_mid   = (r['gr_start'] + r['gr_end']) / 2
            dist = int(abs(atac_mid - gr_mid))
            el['e2g_links'].append({
                'link_type':      link_type,
                'distance_to_tss': dist,
                'correlation':    None,
                'padj':           None,
                'hicar_score':    None,
                'hicar_padj':     None,
            })

    # Strip internal bookkeeping keys and sort by genomic position
    elements = []
    for el in by_peak.values():
        el.pop('_seen_chips', None)
        el.pop('_seen_motifs', None)
        el.pop('_seen_links', None)
        elements.append(el)
    elements.sort(key=lambda e: (e['chr'], e['start']))
    return elements


def _query_tf_binding_peaks(db, tf: str, gene: str) -> list:
    rows = db.execute("""
        SELECT DISTINCT
            tp.chr, tp.start, tp.end, tp.tf_peak_name,
            tp.source, td.dataset, td.cell_type, td.dataset_id
        FROM tf_datasets td
        JOIN tf_peaks          tp  ON td.dataset_id      = tp.dataset_id
        JOIN tf_gene_overlaps  tgo ON tp.peak_id         = tgo.peak_id
        JOIN gene_regions      gr  ON tgo.gene_region_id = gr.gene_region_id
        WHERE td.tf_gene_name = ? AND gr.gene_name = ?
    """, (tf, gene)).fetchall()
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
            tp.chr, tp.start AS tf_start, tp.end AS tf_end,
            ato.overlap_bp
        FROM atac_tf_overlaps ato
        JOIN tf_peaks    tp ON ato.peak_id    = tp.peak_id
        JOIN tf_datasets td ON tp.dataset_id  = td.dataset_id
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
            gr.gene_name, gr.gene_region_subtype, gr.multiome_hicar_support,
            gr.chr AS gr_chr, gr.start AS gr_start, gr.end AS gr_end,
            td.tf_gene_name,
            (ap.start + ap.end) / 2 AS peak_mid,
            (
                SELECT (tss.start + tss.end) / 2
                FROM gene_regions tss
                WHERE tss.gene_name = gr.gene_name
                  AND tss.gene_region_subtype = 'proximal'
                ORDER BY ABS((tss.start + tss.end) / 2 - (ap.start + ap.end) / 2)
                LIMIT 1
            ) AS tss
        FROM atac_tf_overlaps  ato
        JOIN atac_peaks        ap  ON ato.atac_peak_id    = ap.atac_peak_id
        JOIN tf_peaks          tp  ON ato.peak_id         = tp.peak_id
        JOIN tf_datasets       td  ON tp.dataset_id        = td.dataset_id
        JOIN tf_gene_overlaps  tgo ON tp.peak_id           = tgo.peak_id
        JOIN gene_regions      gr  ON tgo.gene_region_id   = gr.gene_region_id
        WHERE ato.atac_peak_id = ?
    """, (atac_peak_id,)).fetchall()
    groups: dict = {}
    for r in rows:
        g = r["gene_name"]
        if g not in groups:
            gr_chr = _norm_chr(str(r["gr_chr"])) if r["gr_chr"] else None
            # tss is the closest proximal-region centre; fall back to linked region midpoint
            tss = r["tss"] if r["tss"] is not None else (r["gr_start"] + r["gr_end"]) // 2
            dist = abs(tss - r["peak_mid"])
            dist_label = f"{dist / 1000:.1f} kb" if dist >= 1000 else f"{dist:,} bp"
            groups[g] = {
                "gene":          g,
                "link_type":     _subtype_to_link_type(r["gene_region_subtype"], r["multiome_hicar_support"]),
                "mediating_tfs": set(),
                "gr_chr":        gr_chr,
                "gr_start":      r["gr_start"],
                "gr_end":        r["gr_end"],
                "tss":           tss,
                "tss_chr":       gr_chr,
                "tss_dist":      dist,
                "tss_dist_label": dist_label,
            }
        groups[g]["mediating_tfs"].add(r["tf_gene_name"])
    result = [dict(g, mediating_tfs=sorted(t for t in g["mediating_tfs"] if t is not None)) for g in groups.values()]
    result.sort(key=lambda g: g["gene"])
    return result


_TP_ORDER = ['ES_0h', 'DE_12h', 'DE_24h', 'DE_36h', 'DE_48h', 'DE_60h', 'DE_72h']

def _query_atac_counts(db, atac_peak_id: int) -> list:
    """Return per-timepoint mean normalized_count + z-score with replicates."""
    try:
        rows = db.execute(
            "SELECT sample, normalized_count, zscore FROM atac_peak_counts WHERE atac_peak_id=?",
            (atac_peak_id,)
        ).fetchall()
    except Exception:
        return []
    by_tp: dict = {}
    for r in rows:
        tp = re.sub(r'_\d+$', '', r['sample'])
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


def _query_gene_elements(db, gene: str) -> list:
    # Two-step to avoid the GROUP BY + COUNT(DISTINCT) over a massive 6-table join.
    # Step 1: get distinct element rows without counting TFs (DISTINCT + no aggregation
    # lets SQLite's hash early-terminate once it has seen each unique atac_peak_id).
    rows = db.execute("""
        SELECT DISTINCT
            ap.atac_peak_id, ap.chr, ap.start, ap.end,
            gr.gene_region_subtype, gr.multiome_hicar_support
        FROM gene_regions      gr
        JOIN tf_gene_overlaps  tgo ON gr.gene_region_id = tgo.gene_region_id
        JOIN tf_peaks          tp  ON tgo.peak_id       = tp.peak_id
        JOIN atac_tf_overlaps  ato ON tp.peak_id        = ato.peak_id
        JOIN atac_peaks        ap  ON ato.atac_peak_id  = ap.atac_peak_id
        WHERE gr.gene_name = ?
        ORDER BY ap.chr, ap.start
    """, (gene,)).fetchall()

    if not rows:
        return []

    # Step 2: for each distinct element, count TFs via the atac_peak_id index.
    # This is a tiny correlated query per element (indexed on atac_peak_id).
    peak_ids = list({r["atac_peak_id"] for r in rows})
    placeholders = ",".join("?" * len(peak_ids))
    tf_counts = {
        r["atac_peak_id"]: r["n_tfs"]
        for r in db.execute(f"""
            SELECT ato.atac_peak_id, COUNT(DISTINCT td.tf_gene_name) AS n_tfs
            FROM atac_tf_overlaps ato
            JOIN tf_peaks         tp  ON ato.peak_id    = tp.peak_id
            JOIN tf_datasets      td  ON tp.dataset_id  = td.dataset_id
            WHERE ato.atac_peak_id IN ({placeholders})
            GROUP BY ato.atac_peak_id
        """, peak_ids).fetchall()
    }

    return [
        {
            "atac_peak_id": r["atac_peak_id"],
            "chr":          _norm_chr(r["chr"]),
            "start":        r["start"],
            "end":          r["end"],
            "link_type":    _subtype_to_link_type(r["gene_region_subtype"], r["multiome_hicar_support"]),
            "n_tfs":        tf_counts.get(r["atac_peak_id"], 0),
        }
        for r in rows
    ]


@perturbseq_bp.route("/api/link-tfs")
def api_link_tfs():
    db = get_db()
    rows = db.execute("SELECT DISTINCT tf_gene_name FROM tf_datasets WHERE tf_gene_name IS NOT NULL ORDER BY tf_gene_name").fetchall()
    return jsonify([r[0] for r in rows])


@perturbseq_bp.route("/api/link-genes")
def api_link_genes():
    db = get_db()
    rows = db.execute("SELECT DISTINCT gene_name FROM gene_regions WHERE gene_name IS NOT NULL ORDER BY gene_name").fetchall()
    return jsonify([r[0] for r in rows])


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
        "SELECT atac_peak_id, chr, start, end FROM atac_peaks WHERE atac_peak_id=?",
        (atac_peak_id,)
    ).fetchone()
    if peak_row is None:
        return None
    peak = dict(peak_row)
    peak["chr"] = _norm_chr(peak["chr"])

    rows = db.execute("""
        SELECT DISTINCT
            td.tf_gene_name,
            td.dataset_id, td.cell_type, td.dataset,
            tp.source AS peak_source,
            gr.gene_region_subtype, gr.multiome_hicar_support,
            gr.start AS gr_start, gr.end AS gr_end,
            ato.overlap_bp
        FROM atac_tf_overlaps  ato
        JOIN tf_peaks          tp  ON ato.peak_id          = tp.peak_id
        JOIN tf_datasets       td  ON tp.dataset_id        = td.dataset_id
        JOIN tf_gene_overlaps  tgo ON tp.peak_id           = tgo.peak_id
        JOIN gene_regions      gr  ON tgo.gene_region_id   = gr.gene_region_id
        WHERE ato.atac_peak_id = ? AND gr.gene_name = ?
        ORDER BY td.tf_gene_name, td.cell_type
    """, (atac_peak_id, gene)).fetchall()

    if not rows:
        return None

    link_type = _subtype_to_link_type(rows[0]["gene_region_subtype"], rows[0]["multiome_hicar_support"])
    peak_mid = (peak["start"] + peak["end"]) // 2
    gr_mid   = (rows[0]["gr_start"] + rows[0]["gr_end"]) // 2
    dist     = abs(peak_mid - gr_mid)
    dist_label = f"{dist / 1000:.1f} kb" if dist >= 1000 else f"{dist:,} bp"

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
        "link_type":    link_type,
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
    gene_expr   = rows_to_dicts(db.execute(
        f"SELECT timepoint, mean_tpm FROM gene_expression WHERE gene=? {_tp_order()}",
        (gene,)))
    return render_template(
        "perturbseq/element_gene_link.html",
        peak=data["peak"],
        gene=gene,
        link_type=data["link_type"],
        dist_label=data["dist_label"],
        mediating_tfs=data["mediating_tfs"],
        atac_counts=atac_counts,
        gene_expr=gene_expr,
    )


@perturbseq_bp.route("/element/<int:atac_peak_id>")
def element_page(atac_peak_id):
    db = get_db()
    row = db.execute(
        "SELECT atac_peak_id, chr, start, end, atac_peak_name FROM atac_peaks WHERE atac_peak_id=?",
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
        """SELECT chr, start, end, atac_peak_name
           FROM atac_peaks
           WHERE chr = ? AND start >= ? AND end <= ? AND atac_peak_id != ?
           ORDER BY start""",
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
    for el in elements:
        el['atac_counts'] = _query_atac_counts(db, el['atac_peak_id'])
    tf_peaks_data = _query_tf_binding_peaks(db, tf, gene)

    tf_expr = rows_to_dicts(db.execute(
        f"SELECT timepoint, mean_tpm FROM gene_expression WHERE gene=? {_tp_order()}",
        (tf,)))
    gene_expr = rows_to_dicts(db.execute(
        f"SELECT timepoint, mean_tpm FROM gene_expression WHERE gene=? {_tp_order()}",
        (gene,)))

    tf_desc_row = db.execute(
        "SELECT name FROM tf_descriptions WHERE tf=?", (tf,)).fetchone()
    tf_desc = tf_desc_row['name'] if tf_desc_row else None

    gene_desc_row = db.execute(
        "SELECT name FROM tf_descriptions WHERE tf=?", (gene,)).fetchone()
    gene_desc = gene_desc_row['name'] if gene_desc_row else None

    return render_template(
        "perturbseq/tf_gene_link.html",
        tf=tf,
        gene=gene,
        tf_desc=tf_desc,
        gene_desc=gene_desc,
        elements=elements,
        tf_peaks_data=tf_peaks_data,
        tf_expr=tf_expr,
        gene_expr=gene_expr,
    )


@perturbseq_bp.route("/dataset/<int:dataset_id>")
def dataset_page(dataset_id):
    db = get_db()

    ds = db.execute(
        "SELECT * FROM tf_datasets WHERE dataset_id=?", (dataset_id,)
    ).fetchone()
    if not ds:
        abort(404)

    tf_desc = db.execute(
        "SELECT name, summary FROM tf_descriptions WHERE tf=?", (ds["tf_gene_name"],)
    ).fetchone()

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
            """SELECT COUNT(DISTINCT gr.gene_name)
               FROM tf_gene_overlaps tgo
               JOIN tf_peaks tp ON tgo.peak_id = tp.peak_id
               JOIN gene_regions gr ON tgo.gene_region_id = gr.gene_region_id
               WHERE tp.dataset_id=?""", (dataset_id,)
        ).fetchone()[0],
    }

    gene_rows = rows_to_dicts(db.execute("""
        SELECT gr.gene_name, gr.gene_region_subtype,
               COUNT(DISTINCT tp.peak_id) AS peak_count
        FROM tf_gene_overlaps tgo
        JOIN tf_peaks tp      ON tgo.peak_id       = tp.peak_id
        JOIN gene_regions gr  ON tgo.gene_region_id = gr.gene_region_id
        WHERE tp.dataset_id=?
        GROUP BY gr.gene_name, gr.gene_region_subtype
        ORDER BY gr.gene_name
    """, (dataset_id,)))

    element_rows = rows_to_dicts(db.execute("""
        SELECT ap.atac_peak_id, ap.chr, ap.start, ap.end,
               MAX(ato.overlap_bp) AS overlap_bp,
               COUNT(DISTINCT tgo.gene_region_id) AS gene_count
        FROM atac_tf_overlaps ato
        JOIN tf_peaks   tp ON ato.peak_id      = tp.peak_id
        JOIN atac_peaks ap ON ato.atac_peak_id = ap.atac_peak_id
        LEFT JOIN tf_gene_overlaps tgo ON tp.peak_id = tgo.peak_id
        WHERE tp.dataset_id=?
        GROUP BY ap.atac_peak_id
        ORDER BY ap.chr, ap.start
    """, (dataset_id,)))
    for r in element_rows:
        r["chr"] = _norm_chr(str(r["chr"]))

    related = rows_to_dicts(db.execute("""
        SELECT dataset_id, dataset, cell_type, cell_type_group, source
        FROM tf_datasets
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
