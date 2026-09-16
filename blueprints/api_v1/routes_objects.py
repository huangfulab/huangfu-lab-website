"""Object endpoints — one entity each.

Every query here is a plain parameterised SELECT. Where a result set is not
inherently bounded it carries an explicit LIMIT, and the accompanying
`*_truncated` flag tells a caller whether they hit it.
"""
from blueprints.perturbseq_bp import (
    _gene_expr_by_timepoint,
    _module_expr_by_timepoint,
    _query_atac_counts,
)

from .catalog import ENUMS
from .core import (
    QUERY_DEADLINE_OBJECT,
    ApiError,
    api_db,
    api_v1_bp,
    arg_enum,
    arg_enum_list,
    arg_int,
    clean_name,
    collection,
    entity,
    html_link,
    page_params,
    path_link,
    query,
    query_one,
    reject_unknown_args,
    scalar,
)
from .names import (
    GSEA_COLLECTION,
    MFUZZ_SOURCE,
    gene_meta,
    module_display_name,
    resolve_gene,
    resolve_go_accession,
    resolve_module,
    tf_dataset_labels,
)

GO_GENE_LIMIT = 500
GO_LIMIT = 500
ELEMENT_LIMIT = 200
PERTURBATION_LIMIT = 500
COEXPRESSION_MAX = 500
MODULE_GENE_LIMIT = 500
ENRICHMENT_LIMIT = 200
REGULATOR_LIMIT = 500
BINDING_DATASET_LIMIT = 2000
PATHWAY_LIMIT = 500
GRNA_GSEA_LIMIT = 500
MODULE_EXPRESSION_MAX_SIZE = 2500

_HOTSPOT_GSEA = GSEA_COLLECTION["hotspot_supermodule"]
_MFUZZ_GSEA = GSEA_COLLECTION[MFUZZ_SOURCE]

# Fisher tests that separated perfectly return an infinite odds ratio. Bare
# Infinity is invalid JSON, and 1e308 would read as a measured value, so these
# become null with a companion flag.
_ODDS = "CASE WHEN {c}.odds_ratio > 1e308 THEN NULL ELSE {c}.odds_ratio END"
_ODDS_INF = "CASE WHEN {c}.odds_ratio > 1e308 THEN 1 ELSE 0 END"


def _is_perturbed_tf(db, gene_name):
    return scalar(
        db,
        "SELECT 1 FROM gsea_tf_table WHERE gene_name = ? "
        "AND gene_set_collection = 'DE-hotspot_modules' LIMIT 1",
        (gene_name,),
    ) is not None


def _merge_module_evidence(perturbation, binding):
    """Combine the perturbation and binding arms into one edge list keyed by
    (collection, module), tagging each with which evidence supports it."""
    merged = {}
    for row in perturbation:
        key = (row["module_collection"], row["module"])
        merged[key] = {
            "module_collection": row["module_collection"],
            "module": module_display_name(row["module"], row["module_collection"]),
            "module_id": None,
            "evidence": "perturbation",
            "direction": row["direction"],
            "mean_NES": row["mean_NES"],
            "n_grnas": row["n_grnas"],
            "perturbation_padj": row["min_padj"],
            "odds_ratio": None,
            "odds_ratio_infinite": False,
            "binding_padj": None,
            "n_overlap": None,
        }
    for row in binding:
        key = (row["module_collection"], row["module_name"])
        existing = merged.get(key)
        if existing is None:
            existing = {
                "module_collection": row["module_collection"],
                "module": module_display_name(row["module_name"], row["module_collection"]),
                "module_id": row["module_id"],
                "evidence": "binding",
                "direction": None,
                "mean_NES": None,
                "n_grnas": None,
                "perturbation_padj": None,
            }
            merged[key] = existing
        else:
            existing["evidence"] = "both"
            existing["module_id"] = row["module_id"]
        existing["odds_ratio"] = row["odds_ratio"]
        existing["odds_ratio_infinite"] = bool(row["odds_ratio_infinite"])
        existing["binding_padj"] = row["padj_fisher"]
        existing["n_overlap"] = row["n_overlap"]
    return sorted(
        merged.values(),
        key=lambda r: (-abs(r["mean_NES"] or 0), -(r["odds_ratio"] or 0), r["module"]),
    )


@api_v1_bp.route("/gene/<gene>")
def gene(gene):
    reject_unknown_args("include", "coexpression_limit")
    include = arg_enum_list("include", ENUMS["gene_include"])
    coexp_limit = arg_int("coexpression_limit", 100, minimum=1, maximum=COEXPRESSION_MAX)

    db = api_db(QUERY_DEADLINE_OBJECT)
    ref = resolve_gene(db, gene)
    gene_id = ref["gene_id"]

    row = query_one(
        db,
        "SELECT g.gene_id, g.ensg_id, g.gene_name, g.gene_biotype, g.chr, "
        "       g.chrom_start, g.chrom_end, g.description, "
        "       g.Other_designations AS other_designations, g.Summary AS summary, "
        "       g.in_perturbation_library, g.primary_gene, pc.publication_count "
        "FROM gene_table g "
        "LEFT JOIN gene_publication_count pc ON pc.gene_id = g.gene_id "
        "WHERE g.gene_id = ?",
        (gene_id,),
    )

    data = dict(row)
    data["id"] = row["gene_name"]
    data["type"] = "gene"
    data["in_perturbation_library"] = bool(row["in_perturbation_library"])
    data["is_primary_symbol"] = bool(row["primary_gene"])
    data.pop("primary_gene", None)
    data["is_tf"] = _is_perturbed_tf(db, row["gene_name"])

    data["synonyms"] = query(
        db,
        "SELECT synonym, synonym_type FROM gene_synonym WHERE gene_id = ? "
        "ORDER BY synonym_type, synonym",
        (gene_id,),
    )
    data["tss"] = query(
        db,
        "SELECT tss_id, chr, tss_position, strand, promoter_source "
        "FROM gene_tss_table WHERE gene_id = ? ORDER BY tss_position",
        (gene_id,),
    )
    data["expression"] = _gene_expr_by_timepoint(db, gene_id)

    modules = query(
        db,
        "SELECT m.module_id, m.module_name, m.source, m.size "
        "FROM gene_module_table gm JOIN module_table m ON m.module_id = gm.module_id "
        "WHERE gm.gene_id = ? ORDER BY m.source, m.module_name",
        (gene_id,),
    )
    for m in modules:
        m["display_name"] = module_display_name(m["module_name"], m["source"])
    data["modules"] = [m for m in modules if m["module_name"] != "unassigned"]
    data["supermodule"] = next(
        (m["module_name"] for m in data["modules"] if m["source"] == "hotspot_supermodule"), None
    )
    data["submodule"] = next(
        (m["module_name"] for m in data["modules"] if m["source"] == "hotspot_submodule"), None
    )
    if data["submodule"] and not data["supermodule"] and "." in data["submodule"]:
        data["supermodule"] = data["submodule"].rsplit(".", 1)[0]

    go_terms = query(
        db,
        "SELECT t.go_accession, t.go_name, t.namespace, gg.evidence, gg.qualifier "
        "FROM go_gene_table gg JOIN go_term_table t ON t.go_id = gg.go_id "
        "WHERE gg.gene_id = ? ORDER BY t.namespace, t.go_name LIMIT ?",
        (gene_id, GO_LIMIT),
    )
    data["go_terms"] = go_terms
    data["go_terms_truncated"] = len(go_terms) == GO_LIMIT

    data["grnas"] = query(
        db,
        "SELECT grna_id, grna_name, active FROM grna_table WHERE gene_id = ? ORDER BY grna_name",
        (gene_id,),
    )

    effects = query(
        db,
        "SELECT gr.gene_name AS perturbed_gene, gr.grna_id, gr.grna_name, gr.active, "
        "       d.coef, d.z_coef, d.padj "
        "FROM de_results d JOIN grna_table gr ON gr.grna_id = d.grna_id "
        "WHERE d.gene_id = ? ORDER BY d.coef LIMIT ?",
        (gene_id, PERTURBATION_LIMIT),
    )
    data["perturbation_effects"] = effects
    data["perturbation_effects_truncated"] = len(effects) == PERTURBATION_LIMIT

    elements = query(
        db,
        "SELECT ap.atac_peak_id, ap.atac_peak_name, ap.chr, ap.chrom_start, ap.chrom_end, "
        "       MIN(l.distance_bp) AS min_distance_bp, "
        "       MAX(l.link_type)   AS link_type, "
        "       MAX(l.cell_type)   AS multiome_cell_type, "
        "       MIN(l.padj)        AS multiome_padj, "
        "       pgc.correlation    AS peak_gene_correlation "
        "FROM ( "
        "    SELECT atac_peak_id, distance_bp, NULL AS link_type, NULL AS cell_type, NULL AS padj "
        "      FROM atac_tss_links WHERE gene_id = ? "
        "    UNION ALL "
        "    SELECT atac_peak_id, distance_to_tss, link_type, cell_type, padj "
        "      FROM multiome_atac_overlaps WHERE gene_id = ? "
        ") l "
        "JOIN atac_peak_table ap ON ap.atac_peak_id = l.atac_peak_id "
        "LEFT JOIN peak_gene_correlation pgc "
        "       ON pgc.atac_peak_id = l.atac_peak_id AND pgc.gene_id = ? "
        "GROUP BY ap.atac_peak_id "
        "ORDER BY ABS(COALESCE(MIN(l.distance_bp), 1000000000)) "
        "LIMIT ?",
        (gene_id, gene_id, gene_id, ELEMENT_LIMIT),
    )
    data["elements"] = elements
    data["elements_truncated"] = len(elements) == ELEMENT_LIMIT

    if "coexpression" in include:
        # The table is asymmetric, so both directions are needed. ORDER BY and
        # LIMIT must stay inside the subquery: sorting after the gene_table join
        # turns a 0.05 s query into a 2 s one.
        data["coexpression"] = query(
            db,
            "SELECT g.gene_id, g.gene_name, x.z_score FROM ( "
            "    SELECT gene_id_2 AS gid, z_score FROM gene_coexpression WHERE gene_id_1 = ? "
            "    UNION ALL "
            "    SELECT gene_id_1,       z_score FROM gene_coexpression WHERE gene_id_2 = ? "
            "    ORDER BY z_score DESC LIMIT ? "
            ") x JOIN gene_table g ON g.gene_id = x.gid "
            "ORDER BY x.z_score DESC",
            (gene_id, gene_id, coexp_limit),
        )

    links = {
        "self": path_link("gene", row["gene_name"]),
        "html": html_link("gene", row["gene_name"]),
    }
    if data["is_tf"]:
        links["tf"] = path_link("tf", row["gene_name"])
    if data["supermodule"]:
        links["supermodule"] = path_link("module", data["supermodule"])

    return entity(data, links, gene_meta(ref))


@api_v1_bp.route("/tf/<tf>")
def tf(tf):
    reject_unknown_args("include")
    include = arg_enum_list("include", ENUMS["tf_include"])

    db = api_db(QUERY_DEADLINE_OBJECT)
    ref = resolve_gene(db, tf)
    gene_id, gene_name = ref["gene_id"], ref["gene_name"]

    is_perturbed = _is_perturbed_tf(db, gene_name)
    labels = tf_dataset_labels(db, gene_id)
    if not is_perturbed and not labels:
        raise_not_a_tf(gene_name)

    row = query_one(
        db,
        "SELECT gene_id, ensg_id, gene_name, description, in_perturbation_library "
        "FROM gene_table WHERE gene_id = ?",
        (gene_id,),
    )
    data = dict(row)
    data["id"] = gene_name
    data["type"] = "tf"
    data["in_perturbation_library"] = bool(row["in_perturbation_library"])
    data["is_perturbed_tf"] = is_perturbed
    data["is_binding_tf"] = bool(labels)
    data["dataset_tf_labels"] = labels

    perturbation = query(
        db,
        "SELECT module_collection, module, gene_set, direction, n_grnas, mean_NES, min_padj "
        "FROM gsea_tf_table "
        "WHERE gene_name = ? AND gene_set_collection IN (?, ?) "
        "  AND module IS NOT NULL AND module <> 'unassigned' "
        "ORDER BY ABS(mean_NES) DESC LIMIT ?",
        (gene_name, _HOTSPOT_GSEA, _MFUZZ_GSEA, REGULATOR_LIMIT),
    )

    binding = []
    if labels:
        placeholders = ",".join("?" * len(labels))
        binding = query(
            db,
            "SELECT m.module_id, m.module_name, m.source AS module_collection, "
            f"       MAX({_ODDS.format(c='e')})     AS odds_ratio, "
            f"       MAX({_ODDS_INF.format(c='e')}) AS odds_ratio_infinite, "
            "       MIN(e.padj_fisher)   AS padj_fisher, "
            "       MAX(e.n_overlap)     AS n_overlap, "
            "       MAX(e.n_tf_targets)  AS n_tf_targets, "
            "       MAX(e.n_module_genes) AS n_module_genes "
            "FROM tf_module_enrichment e "
            "JOIN module_table m ON m.module_id = e.module_id "
            f"WHERE e.tf_gene_name IN ({placeholders}) AND m.module_name <> 'unassigned' "
            "GROUP BY m.module_id ORDER BY odds_ratio DESC LIMIT ?",
            (*labels, REGULATOR_LIMIT),
        )

    merged = _merge_module_evidence(perturbation, binding)
    data["module_regulation"] = merged[:REGULATOR_LIMIT]
    data["module_regulation_truncated"] = (
        len(merged) > REGULATOR_LIMIT
        or len(perturbation) == REGULATOR_LIMIT
        or len(binding) == REGULATOR_LIMIT
    )

    data["grnas"] = query(
        db,
        "SELECT grna_id, grna_name, active FROM grna_table WHERE gene_id = ? ORDER BY grna_name",
        (gene_id,),
    )
    data["n_binding_datasets"] = scalar(
        db, "SELECT COUNT(*) FROM tf_dataset_gene_table WHERE gene_id = ?", (gene_id,), default=0
    )

    if "binding_datasets" in include:
        datasets = query(
            db,
            "SELECT td.dataset_id, td.dataset, td.tf_gene_name AS dataset_tf_label, "
            "       td.cell_type, td.cell_type_group, td.source, td.extra, tdg.match_type "
            "FROM tf_dataset_gene_table tdg "
            "JOIN tf_dataset_table td ON td.dataset_id = tdg.dataset_id "
            "WHERE tdg.gene_id = ? ORDER BY td.source, td.dataset LIMIT ?",
            (gene_id, BINDING_DATASET_LIMIT),
        )
        data["binding_datasets"] = datasets
        data["binding_datasets_truncated"] = len(datasets) == BINDING_DATASET_LIMIT

    if "pathway_enrichment" in include:
        pathways = query(
            db,
            "SELECT gene_set_collection, gene_set, go_id, direction, n_grnas, mean_NES, min_padj "
            "FROM gsea_tf_table WHERE gene_name = ? AND module_collection IS NULL "
            "ORDER BY min_padj LIMIT ?",
            (gene_name, PATHWAY_LIMIT),
        )
        data["pathway_enrichment"] = pathways
        data["pathway_enrichment_truncated"] = len(pathways) == PATHWAY_LIMIT

    links = {
        "self": path_link("tf", gene_name),
        "html": html_link("gene", gene_name),
        "gene": path_link("gene", gene_name),
    }
    return entity(data, links, gene_meta(ref))


def raise_not_a_tf(gene_name):
    from .core import ApiError

    raise ApiError(
        404,
        f"{gene_name} is not a transcription factor in this dataset — it is neither perturbed "
        f"in the screen nor associated with a binding dataset. See /gene/{gene_name}.",
        "tf_not_found",
    )


@api_v1_bp.route("/module/<module>")
def module(module):
    return _module_payload(
        module, known_args=("source", "genes_offset", "genes_limit", "include")
    )


def _module_payload(module, forced_source=None, known_args=(), self_route="module"):
    reject_unknown_args(*known_args)
    source = forced_source or arg_enum("source", ENUMS["module_source"])
    genes_offset = arg_int("genes_offset", 0, minimum=0)
    genes_limit = arg_int("genes_limit", MODULE_GENE_LIMIT, minimum=1, maximum=MODULE_GENE_LIMIT)
    include = arg_enum_list("include", ENUMS["module_include"])

    db = api_db(QUERY_DEADLINE_OBJECT)
    mod = resolve_module(db, module, source=source)
    module_id, module_name, module_source = mod["module_id"], mod["module_name"], mod["source"]

    description = query_one(
        db,
        "SELECT title, standard, extended FROM module_description WHERE module_id = ?",
        (module_id,),
    ) or {}

    data = {
        "id": mod["display_name"],
        "type": "module",
        "module_id": module_id,
        "module_name": module_name,
        "display_name": mod["display_name"],
        "aliases": mod["aliases"],
        "source": module_source,
        "size": mod["size"],
        "title": description.get("title"),
        "description": description.get("standard"),
        "description_extended": description.get("extended"),
    }

    data["n_genes"] = scalar(
        db, "SELECT COUNT(*) FROM gene_module_table WHERE module_id = ?", (module_id,), default=0
    )
    data["genes"] = query(
        db,
        "SELECT g.gene_id, g.gene_name, g.gene_biotype, g.in_perturbation_library, "
        "       pc.publication_count "
        "FROM gene_module_table gm "
        "JOIN gene_table g ON g.gene_id = gm.gene_id "
        "LEFT JOIN gene_publication_count pc ON pc.gene_id = g.gene_id "
        "WHERE gm.module_id = ? ORDER BY g.gene_name, g.gene_id LIMIT ? OFFSET ?",
        (module_id, genes_limit, genes_offset),
    )
    data["genes_offset"] = genes_offset
    data["genes_limit"] = genes_limit

    enrichment = query(
        db,
        "SELECT term_id, source, go_id, term_name, p_value, term_size, query_size, "
        "       intersection_size, term_precision, recall "
        "FROM go_module_enrichment WHERE module_id = ? ORDER BY p_value LIMIT ?",
        (module_id, ENRICHMENT_LIMIT),
    )
    data["enrichment"] = enrichment
    data["enrichment_truncated"] = len(enrichment) == ENRICHMENT_LIMIT

    perturbation = query(
        db,
        "SELECT gene_name AS tf_gene_name, gene_id, direction, n_grnas, mean_NES, min_padj "
        "FROM gsea_tf_table "
        "WHERE gene_set_collection = ? AND module_collection = ? AND module = ? "
        "ORDER BY ABS(mean_NES) DESC LIMIT ?",
        (GSEA_COLLECTION[module_source], module_source, module_name, REGULATOR_LIMIT),
    )
    # Driven off module_id deliberately: adding a gene_set_collection predicate
    # flips the plan to a full scan of this 521K-row table.
    binding = query(
        db,
        "SELECT tf_gene_name, "
        f"       MAX({_ODDS.format(c='tf_module_enrichment')})     AS odds_ratio, "
        f"       MAX({_ODDS_INF.format(c='tf_module_enrichment')}) AS odds_ratio_infinite, "
        "       MIN(padj_fisher)    AS padj_fisher, "
        "       MAX(n_overlap)      AS n_overlap, "
        "       MAX(n_tf_targets)   AS n_tf_targets, "
        "       MAX(n_module_genes) AS n_module_genes, "
        "       MAX(n_total)        AS n_total "
        "FROM tf_module_enrichment WHERE module_id = ? "
        "GROUP BY tf_gene_name ORDER BY odds_ratio DESC LIMIT ?",
        (module_id, REGULATOR_LIMIT),
    )
    merged = _merge_tf_regulators(perturbation, binding)
    data["tf_regulators"] = merged[:REGULATOR_LIMIT]
    data["tf_regulators_truncated"] = (
        len(merged) > REGULATOR_LIMIT
        or len(perturbation) == REGULATOR_LIMIT
        or len(binding) == REGULATOR_LIMIT
    )

    if module_source == "hotspot_supermodule":
        data["child_submodules"] = query(
            db,
            "SELECT module_id, module_name, size FROM module_table "
            "WHERE source = 'hotspot_submodule' AND module_name GLOB ? ORDER BY module_name",
            (f"{module_name}.*",),
        )
        data["parent_module"] = None
    elif module_source == "hotspot_submodule" and "." in module_name:
        data["child_submodules"] = []
        data["parent_module"] = module_name.rsplit(".", 1)[0]
    else:
        data["child_submodules"] = []
        data["parent_module"] = None

    if "grna_gsea" in include:
        rows = query(
            db,
            "SELECT gr.gene_name AS perturbed_gene, gr.grna_id, gr.grna_name, gr.active, "
            "       x.NES, x.pval, x.padj, x.size "
            "FROM gsea_grna_table x JOIN grna_table gr ON gr.grna_id = x.grna_id "
            "WHERE x.module_collection = ? AND x.module = ? ORDER BY x.padj LIMIT ?",
            (module_source, module_name, GRNA_GSEA_LIMIT),
        )
        data["grna_gsea"] = rows
        data["grna_gsea_truncated"] = len(rows) == GRNA_GSEA_LIMIT

    if "expression" in include:
        if data["n_genes"] > MODULE_EXPRESSION_MAX_SIZE:
            data["expression"] = None
            data["expression_note"] = (
                f"Skipped: this module has {data['n_genes']} genes, over the "
                f"{MODULE_EXPRESSION_MAX_SIZE} limit for the aggregate expression profile."
            )
        else:
            data["expression"] = _module_expr_by_timepoint(db, module_id)

    links = {
        "self": path_link(self_route, mod["display_name"]),
        "html": html_link("module", mod["display_name"]),
    }
    if data["parent_module"]:
        links["parent"] = path_link("module", data["parent_module"])
    return entity(data, links)


@api_v1_bp.route("/submodule/<submodule>")
def submodule(submodule):
    # Same payload as /module, but the collection is pinned so the name can never
    # be ambiguous — 'unassigned' exists in two collections.
    return _module_payload(submodule, forced_source="hotspot_submodule",
                           known_args=("genes_offset", "genes_limit", "include"),
                           self_route="submodule")


@api_v1_bp.route("/go-term/<go_term>")
def go_term(go_term):
    reject_unknown_args("genes_offset", "genes_limit", "include")
    genes_offset = arg_int("genes_offset", 0, minimum=0)
    genes_limit = arg_int("genes_limit", 200, minimum=1, maximum=GO_GENE_LIMIT)
    include = arg_enum_list("include", ENUMS["go_include"])

    db = api_db(QUERY_DEADLINE_OBJECT)
    accession = resolve_go_accession(go_term)
    row = query_one(
        db,
        "SELECT go_id, go_accession, go_name, namespace, definition, is_obsolete "
        "FROM go_term_table WHERE go_accession = ?",
        (accession,),
    )
    if not row:
        raise ApiError(404, f"No GO term matching {go_term!r}.", "go_term_not_found")

    go_id = row["go_id"]
    data = dict(row)
    data["id"] = row["go_accession"]
    data["type"] = "go_term"
    data["is_obsolete"] = bool(row["is_obsolete"])
    data["n_genes"] = scalar(
        db, "SELECT COUNT(*) FROM go_gene_table WHERE go_id = ?", (go_id,), default=0
    )
    data["genes"] = query(
        db,
        "SELECT g.gene_id, g.gene_name, g.gene_biotype, gg.evidence, gg.qualifier "
        "FROM go_gene_table gg JOIN gene_table g ON g.gene_id = gg.gene_id "
        "WHERE gg.go_id = ? ORDER BY g.gene_name, g.gene_id LIMIT ? OFFSET ?",
        (go_id, genes_limit, genes_offset),
    )
    data["genes_offset"] = genes_offset
    data["genes_limit"] = genes_limit

    if "module_enrichment" in include:
        data["module_enrichment"] = query(
            db,
            "SELECT e.module_id, m.module_name, m.source, e.term_id, e.source AS term_source, "
            "       e.p_value, e.term_size, e.query_size, e.intersection_size, "
            "       e.term_precision, e.recall "
            "FROM go_module_enrichment e JOIN module_table m ON m.module_id = e.module_id "
            "WHERE e.go_id = ? ORDER BY e.p_value LIMIT ?",
            (go_id, REGULATOR_LIMIT),
        )
    if "tf_enrichment" in include:
        data["tf_enrichment"] = query(
            db,
            "SELECT gene_name AS tf_gene_name, gene_id, gene_set_collection, gene_set, "
            "       direction, n_grnas, mean_NES, min_padj "
            "FROM gsea_tf_table WHERE go_id = ? ORDER BY min_padj LIMIT ?",
            (go_id, REGULATOR_LIMIT),
        )

    links = {
        "self": path_link("go-term", row["go_accession"]),
        "html": html_link("go", row["go_accession"]),
    }
    return entity(data, links)


@api_v1_bp.route("/dataset/<namespace>/<int:dataset_id>")
def dataset(namespace, dataset_id):
    reject_unknown_args()
    if namespace not in ENUMS["dataset_namespace"]:
        raise ApiError(
            404,
            f"Unknown dataset namespace {namespace!r}. Use one of "
            f"{', '.join(ENUMS['dataset_namespace'])}.",
            "dataset_not_found",
        )

    db = api_db(QUERY_DEADLINE_OBJECT)
    if namespace == "tf":
        row = query_one(
            db,
            "SELECT td.dataset_id, td.dataset, td.tf_gene_name AS dataset_tf_label, "
            "       td.cell_type, td.cell_type_group, td.extra, td.source, "
            "       tdg.gene_id AS tf_gene_id, tdg.match_type, g.gene_name AS tf_gene_name "
            "FROM tf_dataset_table td "
            "LEFT JOIN tf_dataset_gene_table tdg ON tdg.dataset_id = td.dataset_id "
            "LEFT JOIN gene_table g ON g.gene_id = tdg.gene_id "
            "WHERE td.dataset_id = ?",
            (dataset_id,),
        )
        if not row:
            raise ApiError(404, f"No TF dataset with id {dataset_id}.", "dataset_not_found")
        data = dict(row)
        data["id"] = dataset_id
        data["type"] = "dataset"
        data["namespace"] = "tf"
        # Covering index on tf_peaks(dataset_id); the peak list itself is never
        # served — ~6,900 rows per dataset over a 121M-row table.
        data["n_tf_peaks"] = scalar(
            db, "SELECT COUNT(*) FROM tf_peaks WHERE dataset_id = ?", (dataset_id,), default=0
        )
        links = {"self": path_link("dataset", "tf", dataset_id)}
        if row["tf_gene_name"]:
            links["tf"] = path_link("tf", row["tf_gene_name"])
        return entity(data, links)

    row = query_one(
        db,
        "SELECT ptm_dataset_id, dataset, source, ptm_name, cell_type, timepoint "
        "FROM ptm_dataset_table WHERE ptm_dataset_id = ?",
        (dataset_id,),
    )
    if not row:
        raise ApiError(404, f"No PTM dataset with id {dataset_id}.", "dataset_not_found")
    data = dict(row)
    data["id"] = dataset_id
    data["type"] = "dataset"
    data["namespace"] = "ptm"
    data["n_atac_peaks"] = scalar(
        db, "SELECT COUNT(*) FROM ptm_atac_links WHERE ptm_dataset_id = ?",
        (dataset_id,), default=0,
    )
    return entity(data, {"self": path_link("dataset", "ptm", dataset_id)})


def _resolve_peak(db, peak):
    if peak.isdigit():
        row = query_one(
            db,
            "SELECT atac_peak_id, chr, chrom_start, chrom_end, atac_peak_name "
            "FROM atac_peak_table WHERE atac_peak_id = ?",
            (int(peak),),
        )
    else:
        row = query_one(
            db,
            "SELECT atac_peak_id, chr, chrom_start, chrom_end, atac_peak_name "
            "FROM atac_peak_table WHERE atac_peak_name = ? LIMIT 1",
            (clean_name(peak, "peak"),),
        )
    if not row:
        raise ApiError(404, f"No ATAC peak matching {peak!r}.", "peak_not_found")
    return row


@api_v1_bp.route("/atac-peak/<peak>")
def atac_peak(peak):
    reject_unknown_args()
    db = api_db(QUERY_DEADLINE_OBJECT)
    row = _resolve_peak(db, peak)
    peak_id = row["atac_peak_id"]

    data = dict(row)
    data["id"] = peak_id
    data["type"] = "atac_peak"
    data["width"] = row["chrom_end"] - row["chrom_start"]

    linked = query(
        db,
        "SELECT g.gene_id, g.gene_name, l.distance_bp, l.link_type, l.cell_type, "
        "       l.beta, l.z, l.p, l.padj, pgc.correlation "
        "FROM ( "
        "    SELECT gene_id, distance_bp, NULL AS link_type, NULL AS cell_type, "
        "           NULL AS beta, NULL AS z, NULL AS p, NULL AS padj "
        "      FROM atac_tss_links WHERE atac_peak_id = ? "
        "    UNION ALL "
        "    SELECT gene_id, distance_to_tss, link_type, cell_type, beta, z, p, padj "
        "      FROM multiome_atac_overlaps WHERE atac_peak_id = ? "
        ") l "
        "JOIN gene_table g ON g.gene_id = l.gene_id "
        "LEFT JOIN peak_gene_correlation pgc "
        "       ON pgc.atac_peak_id = ? AND pgc.gene_id = l.gene_id "
        "ORDER BY ABS(COALESCE(l.distance_bp, 1000000000)) LIMIT ?",
        (peak_id, peak_id, peak_id, ELEMENT_LIMIT),
    )
    data["linked_genes"] = linked
    data["linked_genes_truncated"] = len(linked) == ELEMENT_LIMIT
    data["n_linked_genes"] = len(linked)

    data["accessibility"] = _query_atac_counts(db, peak_id)
    data["histone_ptms"] = query(
        db,
        "SELECT p.ptm_dataset_id, p.ptm_name, p.cell_type, p.timepoint, p.dataset, p.source "
        "FROM ptm_atac_links l JOIN ptm_dataset_table p ON p.ptm_dataset_id = l.ptm_dataset_id "
        "WHERE l.atac_peak_id = ? ORDER BY p.ptm_name",
        (peak_id,),
    )
    # Covering index, ~0.3 ms. The TF list itself is a separate endpoint: a busy
    # peak overlaps 10,000+ TF peaks and GROUP BY defeats any LIMIT.
    data["n_tf_peak_overlaps"] = scalar(
        db, "SELECT COUNT(*) FROM atac_tf_overlaps WHERE atac_peak_id = ?", (peak_id,), default=0
    )

    links = {
        "self": path_link("atac-peak", peak_id),
        "tfs": path_link("atac-peak", peak_id) + "/tfs",
    }
    return entity(data, links)


@api_v1_bp.route("/atac-peak/<peak>/tfs")
def atac_peak_tfs(peak):
    reject_unknown_args("page", "per_page")
    page, per_page, offset = page_params()

    db = api_db(QUERY_DEADLINE_OBJECT)
    peak_id = _resolve_peak(db, peak)["atac_peak_id"]
    # CROSS JOIN pins entry at atac_tf_overlaps' peak index, then walks rowids.
    # Any other order scans a 122M-row table.
    rows = query(
        db,
        "SELECT td.tf_gene_name, g.gene_id AS tf_gene_id, td.source AS dataset_source, "
        "       COUNT(DISTINCT td.dataset_id) AS n_datasets, "
        "       MAX(ato.overlap_bp) AS max_overlap_bp "
        "FROM atac_tf_overlaps ato "
        "CROSS JOIN tf_peaks tp ON tp.peak_id = ato.peak_id "
        "CROSS JOIN tf_dataset_table td ON td.dataset_id = tp.dataset_id "
        "LEFT JOIN tf_dataset_gene_table tdg ON tdg.dataset_id = td.dataset_id "
        "LEFT JOIN gene_table g ON g.gene_id = tdg.gene_id "
        "WHERE ato.atac_peak_id = ? "
        "GROUP BY td.tf_gene_name, td.source "
        # Orders on the full GROUP BY key, so the sort is total even though no
        # single column here is unique. /* stable-order */ tells the pagination
        # lint this was checked.
        "ORDER BY n_datasets DESC, td.tf_gene_name, td.source /* stable-order */ "
        "LIMIT ? OFFSET ?",
        (peak_id, per_page, offset),
    )
    for r in rows:
        if r["tf_gene_id"]:
            r["url"] = path_link("tf", r["tf_gene_name"])
    # total is deliberately null: counting distinct TFs costs as much as the
    # query itself on a busy peak.
    return collection(rows, page=page, per_page=per_page, total=None, total_is_exact=False)


def _merge_tf_regulators(perturbation, binding):
    merged = {}
    for row in perturbation:
        merged[row["tf_gene_name"]] = {
            "tf_gene_name": row["tf_gene_name"],
            "gene_id": row["gene_id"],
            "evidence": "perturbation",
            "direction": row["direction"],
            "mean_NES": row["mean_NES"],
            "n_grnas": row["n_grnas"],
            "perturbation_padj": row["min_padj"],
            "odds_ratio": None,
            "odds_ratio_infinite": False,
            "binding_padj": None,
            "n_overlap": None,
        }
    for row in binding:
        existing = merged.get(row["tf_gene_name"])
        if existing is None:
            existing = {
                "tf_gene_name": row["tf_gene_name"],
                "gene_id": None,
                "evidence": "binding",
                "direction": None,
                "mean_NES": None,
                "n_grnas": None,
                "perturbation_padj": None,
            }
            merged[row["tf_gene_name"]] = existing
        else:
            existing["evidence"] = "both"
        existing["odds_ratio"] = row["odds_ratio"]
        existing["odds_ratio_infinite"] = bool(row["odds_ratio_infinite"])
        existing["binding_padj"] = row["padj_fisher"]
        existing["n_overlap"] = row["n_overlap"]
    return sorted(
        merged.values(),
        key=lambda r: (-abs(r["mean_NES"] or 0), -(r["odds_ratio"] or 0), r["tf_gene_name"]),
    )
