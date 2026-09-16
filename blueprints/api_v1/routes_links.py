"""Link endpoints — evidence connecting two entities.

Several queries here use CROSS JOIN to pin the join order. That is not stylistic:
ANALYZE has never been run on this database, so SQLite plans these joins on
default cardinality guesses and picks catastrophically wrong orders on the
100M+ row tables. Each CROSS JOIN below is marked with what it prevents.
"""
from .catalog import ENUMS
from .core import (
    QUERY_DEADLINE_COLLECTION,
    QUERY_DEADLINE_OBJECT,
    ApiError,
    api_db,
    api_v1_bp,
    arg_bool,
    arg_enum,
    arg_enum_list,
    arg_float,
    arg_str,
    collection,
    entity,
    html_link,
    page_params,
    path_link,
    query,
    reject_unknown_args,
    scalar,
)
from .names import (
    GSEA_COLLECTION,
    HIDDEN_MODULES,
    module_display_name,
    resolve_gene,
    resolve_module,
    tf_dataset_labels,
)

BINDING_SCORE_LIMIT = 500
ELEMENT_EVIDENCE_LIMIT = 500
GRNA_LIMIT = 500
EDGE_ARM_LIMIT = 20000

# NULL-guarded filters rather than a WHERE clause assembled at request time, so
# no user value can reach the SQL text. gene_module_table is only 26K rows, so
# giving up index selectivity on these predicates costs nothing measurable.
_GENE_MODULE_FILTERS = (
    "WHERE (? IS NULL OR gm.gene_id = ?) "
    "  AND (? IS NULL OR gm.module_id = ?) "
    "  AND (? IS NULL OR m.source = ?) "
    "  AND (? = 1 OR m.module_name <> 'unassigned') "
)
_GENE_MODULE_COUNT_SQL = (
    "SELECT COUNT(*) FROM gene_module_table gm "
    "JOIN module_table m ON m.module_id = gm.module_id " + _GENE_MODULE_FILTERS
)
_GENE_MODULE_ROWS_SQL = (
    "SELECT g.gene_id, g.gene_name, g.gene_biotype, g.in_perturbation_library, "
    "       m.module_id, m.module_name, m.source, m.size AS module_size "
    "FROM gene_module_table gm "
    "JOIN module_table m ON m.module_id = gm.module_id "
    "JOIN gene_table g ON g.gene_id = gm.gene_id " + _GENE_MODULE_FILTERS +
    "ORDER BY m.module_name, g.gene_name, g.gene_id LIMIT ? OFFSET ?"
)


@api_v1_bp.route("/link/tf-gene/<tf>/<gene>")
def link_tf_gene(tf, gene):
    reject_unknown_args("include")
    include = arg_enum_list("include", ENUMS["tf_gene_include"])

    db = api_db(QUERY_DEADLINE_OBJECT)
    tf_ref = resolve_gene(db, tf)
    gene_ref = resolve_gene(db, gene)
    tf_id, gene_id = tf_ref["gene_id"], gene_ref["gene_id"]

    # Tier 1. The two CROSS JOINs force entry through the (dataset_id, gene_id)
    # primary key of tf_gene_binding_scores. Left to itself the planner picks
    # idx_tf_bs_gene_id and reads ~6,600 rows of a 184M-row table: 4.95 s versus
    # 0.008 s. Driving from tf_dataset_gene_table also bounds the outer loop at
    # the TF's dataset count and resolves aliased TF labels for free.
    datasets = query(
        db,
        "SELECT td.dataset_id, td.dataset, td.tf_gene_name AS dataset_tf_label, "
        "       td.cell_type, td.cell_type_group, td.source, tdg.match_type, "
        "       bs.binding_score_A, bs.binding_score_B, bs.binding_score_C, bs.binding_score_D "
        "FROM tf_dataset_gene_table tdg "
        "CROSS JOIN tf_dataset_table td ON td.dataset_id = tdg.dataset_id "
        "CROSS JOIN tf_gene_binding_scores bs "
        "        ON bs.dataset_id = td.dataset_id AND bs.gene_id = ? "
        "WHERE tdg.gene_id = ? "
        "ORDER BY bs.binding_score_A DESC LIMIT ?",
        (gene_id, tf_id, BINDING_SCORE_LIMIT),
    )

    perturbation = query(
        db,
        "SELECT gr.grna_id, gr.grna_name, gr.active, d.coef, d.z_coef, d.padj "
        "FROM de_results d JOIN grna_table gr ON gr.grna_id = d.grna_id "
        "WHERE d.gene_id = ? AND gr.gene_id = ? ORDER BY d.coef",
        (gene_id, tf_id),
    )

    scores = [d["binding_score_A"] for d in datasets if d["binding_score_A"] is not None]
    data = {
        "id": f"{tf_ref['gene_name']}->{gene_ref['gene_name']}",
        "type": "tf_gene_link",
        "tf_gene_name": tf_ref["gene_name"], "tf_gene_id": tf_id,
        "gene_name": gene_ref["gene_name"], "gene_id": gene_id,
        "binding_evidence": bool(datasets),
        "n_datasets_with_binding": len(datasets),
        "max_binding_score_A": max(scores) if scores else None,
        "datasets": datasets,
        "datasets_truncated": len(datasets) == BINDING_SCORE_LIMIT,
        "perturbation": perturbation,
        "perturbation_evidence": bool(perturbation),
    }

    if "elements" in include:
        if not datasets:
            # The gate. Most of the 1,705 x 24,960 TF-gene space has no binding at
            # all, and the peak-level query costs ~0.7 s even when it will return
            # nothing. Skipping it here is what keeps this endpoint usable.
            data["elements"] = []
            data["elements_note"] = (
                "Skipped: no dataset places this TF at this gene, so there is no peak-level "
                "evidence to enumerate."
            )
        else:
            data["elements"] = query(
                db,
                "WITH gene_peaks AS ( "
                "    SELECT atac_peak_id FROM atac_tss_links WHERE gene_id = ? "
                "    UNION "
                "    SELECT atac_peak_id FROM multiome_atac_overlaps WHERE gene_id = ? "
                ") "
                "SELECT DISTINCT gp.atac_peak_id, ap.atac_peak_name, ap.chr, "
                "       ap.chrom_start, ap.chrom_end, ato.overlap_bp, "
                "       tp.peak_id, tp.source AS peak_source, "
                "       td.dataset_id, td.dataset, td.cell_type "
                "FROM gene_peaks gp "
                # Anchors the walk on the gene's own peaks. Inverted, this scans
                # 122M rows of atac_tf_overlaps.
                "CROSS JOIN atac_tf_overlaps ato ON ato.atac_peak_id = gp.atac_peak_id "
                "CROSS JOIN tf_peaks tp ON tp.peak_id = ato.peak_id "
                "CROSS JOIN tf_dataset_table td ON td.dataset_id = tp.dataset_id "
                "CROSS JOIN atac_peak_table ap ON ap.atac_peak_id = gp.atac_peak_id "
                "WHERE tp.dataset_id IN "
                "      (SELECT dataset_id FROM tf_dataset_gene_table WHERE gene_id = ?) "
                "LIMIT ?",
                (gene_id, gene_id, tf_id, ELEMENT_EVIDENCE_LIMIT),
            )
            data["elements_truncated"] = len(data["elements"]) == ELEMENT_EVIDENCE_LIMIT

    links = {
        "self": path_link("link", "tf-gene", tf_ref["gene_name"], gene_ref["gene_name"]),
        "tf": path_link("tf", tf_ref["gene_name"]),
        "gene": path_link("gene", gene_ref["gene_name"]),
        "html": html_link("link", tf_ref["gene_name"], gene_ref["gene_name"]),
    }
    return entity(data, links)


@api_v1_bp.route("/link/tf-module/<tf>/<module>")
def link_tf_module(tf, module):
    reject_unknown_args("source")
    source = arg_enum("source", ENUMS["module_source"])

    db = api_db(QUERY_DEADLINE_OBJECT)
    tf_ref = resolve_gene(db, tf)
    mod = resolve_module(db, module, source=source)
    tf_id, tf_name = tf_ref["gene_id"], tf_ref["gene_name"]
    module_id, module_name, module_source = mod["module_id"], mod["module_name"], mod["source"]

    perturbation = query(
        db,
        "SELECT direction, n_grnas, mean_NES, min_padj, gene_set "
        "FROM gsea_tf_table WHERE gene_name = ? AND gene_set_collection = ? "
        "  AND module_collection = ? AND module = ?",
        (tf_name, GSEA_COLLECTION[module_source], module_source, module_name),
    )

    labels = tf_dataset_labels(db, tf_id) or [tf_name]
    placeholders = ",".join("?" * len(labels))
    # Anchored on module_id, which is indexed. Filtering this table by
    # gene_set_collection instead would scan all 521K rows, and that column uses
    # 'mfuzz' where module_table.source says 'mfuzz_k7'.
    binding = query(
        db,
        "SELECT tf_gene_name, "
        "       CASE WHEN odds_ratio > 1e308 THEN NULL ELSE odds_ratio END AS odds_ratio, "
        "       CASE WHEN odds_ratio > 1e308 THEN 1 ELSE 0 END AS odds_ratio_infinite, "
        "       pval_fisher, padj_fisher, chisq_stat, pval_chisq, "
        "       n_overlap, n_tf_targets, n_module_genes, n_total "
        "FROM tf_module_enrichment "
        f"WHERE module_id = ? AND tf_gene_name IN ({placeholders}) "
        "ORDER BY padj_fisher LIMIT 1",
        (module_id, *labels),
    )

    per_grna = query(
        db,
        "SELECT gr.grna_id, gr.grna_name, gr.active, x.NES, x.pval, x.padj, x.size "
        "FROM gsea_grna_table x JOIN grna_table gr ON gr.grna_id = x.grna_id "
        "WHERE x.module_collection = ? AND x.module = ? AND gr.gene_id = ? "
        "ORDER BY x.padj LIMIT ?",
        (module_source, module_name, tf_id, GRNA_LIMIT),
    )

    overlap = query(
        db,
        "SELECT COUNT(*) AS n_module_genes_tested, AVG(d.coef) AS mean_coef, "
        "       SUM(CASE WHEN d.padj < 0.05 THEN 1 ELSE 0 END) AS n_significant "
        "FROM gene_module_table gm JOIN de_results d ON d.gene_id = gm.gene_id "
        "WHERE gm.module_id = ? "
        "  AND d.grna_id IN (SELECT grna_id FROM grna_table WHERE gene_id = ?)",
        (module_id, tf_id),
    )

    p = perturbation[0] if perturbation else {}
    b = binding[0] if binding else {}
    o = overlap[0] if overlap else {}
    evidence = (
        "both" if perturbation and binding
        else "perturbation" if perturbation
        else "binding" if binding
        else "none"
    )

    data = {
        "id": f"{tf_name}->{mod['display_name']}",
        "type": "tf_module_link",
        "tf_gene_name": tf_name, "tf_gene_id": tf_id,
        "module_id": module_id, "module_name": module_name,
        "module_display_name": mod["display_name"], "module_collection": module_source,
        "evidence": evidence,
        "direction": p.get("direction"), "mean_NES": p.get("mean_NES"),
        "n_grnas": p.get("n_grnas"), "perturbation_padj": p.get("min_padj"),
        "odds_ratio": b.get("odds_ratio"),
        "odds_ratio_infinite": bool(b.get("odds_ratio_infinite")),
        "binding_padj": b.get("padj_fisher"), "binding_pval": b.get("pval_fisher"),
        "n_overlap": b.get("n_overlap"), "n_tf_targets": b.get("n_tf_targets"),
        "n_module_genes": b.get("n_module_genes"), "n_total": b.get("n_total"),
        "n_module_genes_tested": o.get("n_module_genes_tested"),
        "mean_coef": o.get("mean_coef"), "n_significant": o.get("n_significant"),
        "per_grna": per_grna,
        "per_grna_truncated": len(per_grna) == GRNA_LIMIT,
    }
    links = {
        "self": path_link("link", "tf-module", tf_name, mod["display_name"]),
        "tf": path_link("tf", tf_name),
        "module": path_link("module", mod["display_name"]),
    }
    return entity(data, links)


@api_v1_bp.route("/link/gene-module")
def link_gene_module():
    reject_unknown_args("gene", "module", "source", "include_unassigned", "page", "per_page")
    gene = arg_str("gene")
    module = arg_str("module")
    source = arg_enum("source", ENUMS["module_source"])
    include_unassigned = arg_bool("include_unassigned", False)
    page, per_page, offset = page_params()

    if not (gene or module or source):
        raise ApiError(
            400,
            "Supply at least one of 'gene', 'module' or 'source'. An unfiltered listing would "
            "return every gene-module edge in the dataset.",
            "filter_required",
        )

    db = api_db(QUERY_DEADLINE_COLLECTION)
    gene_id = resolve_gene(db, gene)["gene_id"] if gene else None
    module_id = resolve_module(db, module, source=source)["module_id"] if module else None
    source_filter = source if (source and not module) else None
    params = (
        gene_id, gene_id,
        module_id, module_id,
        source_filter, source_filter,
        1 if include_unassigned else 0,
    )

    total = scalar(
        db, _GENE_MODULE_COUNT_SQL, params, default=0
    )
    rows = query(
        db, _GENE_MODULE_ROWS_SQL, (*params, per_page, offset)
    )
    for r in rows:
        r["module_display_name"] = module_display_name(r["module_name"], r["source"])
        r["in_perturbation_library"] = bool(r["in_perturbation_library"])
    return collection(rows, page=page, per_page=per_page, total=total)


@api_v1_bp.route("/edges")
def edges():
    reject_unknown_args("level", "evidence", "min_odds_ratio", "max_padj", "min_abs_nes",
                        "page", "per_page")
    level = arg_enum("level", ENUMS["module_source"], default="hotspot_supermodule")
    evidence = arg_enum("evidence", ENUMS["edge_evidence"], default="any")
    min_or = arg_float("min_odds_ratio", 1.0, minimum=0)
    max_padj = arg_float("max_padj", 0.05, minimum=0, maximum=1)
    min_nes = arg_float("min_abs_nes", 0.0, minimum=0)
    page, per_page, offset = page_params()

    db = api_db(QUERY_DEADLINE_COLLECTION)

    perturbation = query(
        db,
        "SELECT gene_name AS tf, gene_id AS tf_gene_id, module, module_collection, "
        "       direction, mean_NES, n_grnas, min_padj "
        "FROM gsea_tf_table "
        "WHERE gene_set_collection = ? AND module_collection = ? "
        "  AND module <> 'unassigned' AND ABS(mean_NES) >= ? "
        "ORDER BY ABS(mean_NES) DESC LIMIT ?",
        (GSEA_COLLECTION[level], level, min_nes, EDGE_ARM_LIMIT),
    ) if evidence != "binding" else []

    # Filtered through module_table.source, not tf_module_enrichment's own
    # gene_set_collection column: the latter is unindexed here and turns a
    # 13-module seed into a full scan of 521K rows.
    binding = query(
        db,
        "SELECT e.tf_gene_name AS tf, m.module_name AS module, m.source AS module_collection, "
        "       m.module_id, "
        "       MAX(CASE WHEN e.odds_ratio > 1e308 THEN NULL ELSE e.odds_ratio END) AS odds_ratio, "
        "       MAX(CASE WHEN e.odds_ratio > 1e308 THEN 1 ELSE 0 END) AS odds_ratio_infinite, "
        "       MIN(e.padj_fisher) AS binding_padj, MAX(e.n_overlap) AS n_overlap "
        "FROM module_table m "
        "JOIN tf_module_enrichment e ON e.module_id = m.module_id "
        "WHERE m.source = ? AND m.module_name <> 'unassigned' "
        "  AND e.odds_ratio >= ? AND e.padj_fisher <= ? "
        "GROUP BY e.tf_gene_name, m.module_name, m.source, m.module_id "
        "ORDER BY odds_ratio DESC LIMIT ?",
        (level, min_or, max_padj, EDGE_ARM_LIMIT),
    ) if evidence != "perturbation" else []

    merged = {}
    for r in perturbation:
        if (level, r["module"]) in HIDDEN_MODULES:
            continue
        merged[(r["tf"], r["module"])] = {
            "tf": r["tf"], "tf_gene_id": r["tf_gene_id"],
            "module": module_display_name(r["module"], level),
            "module_id": None, "module_collection": level,
            "evidence": "perturbation", "direction": r["direction"],
            "mean_NES": r["mean_NES"], "n_grnas": r["n_grnas"],
            "perturbation_padj": r["min_padj"],
            "odds_ratio": None, "odds_ratio_infinite": False,
            "binding_padj": None, "n_overlap": None,
        }
    for r in binding:
        if (level, r["module"]) in HIDDEN_MODULES:
            continue
        key = (r["tf"], r["module"])
        entry = merged.get(key)
        if entry is None:
            entry = {
                "tf": r["tf"], "tf_gene_id": None,
                "module": module_display_name(r["module"], level),
                "module_id": r["module_id"], "module_collection": level,
                "evidence": "binding", "direction": None, "mean_NES": None,
                "n_grnas": None, "perturbation_padj": None,
            }
            merged[key] = entry
        else:
            entry["evidence"] = "both"
            entry["module_id"] = r["module_id"]
        entry["odds_ratio"] = r["odds_ratio"]
        entry["odds_ratio_infinite"] = bool(r["odds_ratio_infinite"])
        entry["binding_padj"] = r["binding_padj"]
        entry["n_overlap"] = r["n_overlap"]

    rows = list(merged.values())
    if evidence == "both":
        rows = [r for r in rows if r["evidence"] == "both"]
    elif evidence in ("perturbation", "binding"):
        rows = [r for r in rows if r["evidence"] in (evidence, "both")]

    for r in rows:
        r["strength"] = max(abs(r["mean_NES"] or 0), r["odds_ratio"] or 0)
    rows.sort(key=lambda r: (-r["strength"], r["tf"], r["module"]))

    return collection(
        rows[offset:offset + per_page], page=page, per_page=per_page, total=len(rows),
        meta={"level": level, "evidence_filter": evidence},
    )
