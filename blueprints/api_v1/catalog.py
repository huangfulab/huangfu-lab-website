"""Endpoint catalog for the public API.

Single source of truth, consumed by three places:
  - meta.py            renders the HTML docs page and the JSON index from it
  - core.py            reads ENUMS so validated values and documented values cannot diverge
  - scripts/test_api_v1.py  generates most of its cases from it

Pure data. No imports, no Flask, no DB.
"""

ENUMS = {
    "module_source": ("hotspot_supermodule", "hotspot_submodule", "mfuzz_k7"),
    "regulation_level": ("hotspot_supermodule", "hotspot_submodule", "mfuzz_k7", "all"),
    "gene_include": ("go_terms", "perturbation_effects", "elements", "coexpression"),
    "tf_include": ("binding_datasets", "pathway_enrichment"),
    "module_include": ("genes", "enrichment", "tf_regulators", "grna_gsea", "expression"),
    "tf_set": ("perturbed", "binding", "all"),
    "dataset_namespace": ("tf", "ptm"),
    "edge_evidence": ("any", "perturbation", "binding", "both"),
    "go_include": ("genes", "module_enrichment", "tf_enrichment"),
    "tf_gene_include": ("datasets", "elements"),
    "search_type": ("gene", "tf", "module", "submodule", "gene_cluster", "go_term", "synonym"),
}

# ── reusable param descriptors ───────────────────────────────────────────────

P_PAGE = {
    "name": "page", "type": "integer", "default": 1, "min": 1,
    "description": "1-based page number.",
}
P_PER_PAGE = {
    "name": "per_page", "type": "integer", "default": 25, "min": 1, "max": 500,
    "description": "Rows per page. Values above 500 are rejected, not clamped.",
}
P_MODULE_SOURCE = {
    "name": "source", "type": "enum", "enum": "module_source", "required": False,
    "description": "Module collection. Required when the name is ambiguous — 'unassigned' "
                   "exists in both hotspot_supermodule and hotspot_submodule.",
}

ENDPOINTS = [
    {
        "id": "index",
        "group": "meta",
        "phase": 1,
        "method": "GET",
        "rule": "/",
        "endpoint": "perturbseq_api_v1.api_index",
        "path": "/",
        "summary": "Machine-readable index of every endpoint, with limits and enum values.",
        "path_params": [],
        "query_params": [
            {"name": "examples", "type": "boolean", "default": True,
             "description": "Set false for a compact listing without example response bodies."},
        ],
        "returns": "index",
        "example_url": "/endoderm-perturbseq/api/v1/",
        "errors": [],
        "notes": [],
    },
    {
        "id": "gene",
        "group": "objects",
        "phase": 1,
        "method": "GET",
        "rule": "/gene/<gene>",
        "endpoint": "perturbseq_api_v1.gene",
        "path": "/gene/{gene}",
        "summary": "One gene: identifiers, coordinates, expression over the differentiation "
                   "time course, module membership and gRNAs, plus counts of its GO "
                   "annotations, perturbation effects and linked ATAC peaks. The lists "
                   "themselves are opt-in via include=.",
        "path_params": [
            {"name": "gene", "type": "string", "required": True, "example": "SOX17",
             "description": "HGNC symbol, Ensembl gene ID, or a synonym. Symbols with duplicate "
                            "entries resolve to the primary_gene=1 row."},
        ],
        "query_params": [
            {"name": "include", "type": "enum_list", "enum": "gene_include", "default": "",
             "description": "Comma-separated optional lists, each counted in the default "
                            "response as n_<name>. perturbation_effects is the effect of every "
                            "gRNA in the screen on this gene (~400 rows); coexpression is the "
                            "slowest (~2 s cold)."},
            {"name": "coexpression_limit", "type": "integer", "default": 100, "min": 1, "max": 500,
             "description": "Top co-expressed partners by |z|. Only used with include=coexpression."},
        ],
        "returns": "entity",
        "links": ["self", "html", "tf", "modules"],
        "example_url": "/endoderm-perturbseq/api/v1/gene/SOX17",
        "errors": [
            {"status": 404, "code": "gene_not_found", "when": "No symbol, Ensembl ID or synonym matches."},
        ],
        "notes": [
            "When the input was not the canonical symbol, meta.resolved_from and meta.match_type "
            "report what was matched.",
        ],
    },
    {
        "id": "tf",
        "group": "objects",
        "phase": 1,
        "method": "GET",
        "rule": "/tf/<tf>",
        "endpoint": "perturbseq_api_v1.tf",
        "path": "/tf/{tf}",
        "summary": "One transcription factor: which modules it regulates, by perturbation "
                   "(GSEA over the CRISPR screen) and by binding (ChIP/ATAC enrichment).",
        "path_params": [
            {"name": "tf", "type": "string", "required": True, "example": "ARID1A",
             "description": "Gene symbol of a TF. There is no separate TF table — a TF is a gene "
                            "that is either perturbed in the screen or has binding datasets."},
        ],
        "query_params": [
            {"name": "level", "type": "enum", "enum": "regulation_level",
             "default": "hotspot_supermodule",
             "description": "Which module collection module_regulation covers. The default "
                            "(13 supermodules) keeps the response small; n_module_regulation "
                            "gives the row count for every collection. 'all' returns ~300 rows."},
            {"name": "include", "type": "enum_list", "enum": "tf_include", "default": "",
             "description": "Comma-separated optional blocks. 'binding_datasets' can be large "
                            "(up to ~1,073 datasets for CTCF); 'pathway_enrichment' adds GO/KEGG GSEA."},
        ],
        "returns": "entity",
        "links": ["self", "html", "gene"],
        "example_url": "/endoderm-perturbseq/api/v1/tf/ARID1A",
        "errors": [
            {"status": 404, "code": "tf_not_found",
             "when": "The gene exists but is neither perturbed nor has binding data — the error "
                     "points at /gene/{name}."},
            {"status": 404, "code": "gene_not_found", "when": "No such gene at all."},
        ],
        "notes": [
            "is_perturbed_tf and is_binding_tf are reported separately; they are different "
            "populations (75 perturbed, ~1,705 with binding data).",
            "odds_ratio is null with odds_ratio_infinite=true where the underlying Fisher test "
            "returned an infinite estimate.",
        ],
    },
    {
        "id": "module",
        "group": "objects",
        "phase": 1,
        "method": "GET",
        "rule": "/module/<module>",
        "endpoint": "perturbseq_api_v1.module",
        "path": "/module/{module}",
        "summary": "One co-expression module: title, description, size and hierarchy, plus "
                   "counts of its member genes, enrichment terms and TF regulators. The lists "
                   "themselves are opt-in via include=.",
        "path_params": [
            {"name": "module", "type": "string", "required": True, "example": "DE-1",
             "description": "Supermodule (DE-1), submodule (DE-1.2), or gene cluster (GC1). "
                            "GC1 / TC-1 / gene_cluster_1 / cluster_1 all resolve to the same module."},
        ],
        "query_params": [
            P_MODULE_SOURCE,
            {"name": "include", "type": "enum_list", "enum": "module_include", "default": "",
             "description": "Comma-separated optional lists. 'genes' is paged by genes_offset and "
                            "genes_limit; 'tf_regulators' merges perturbation and binding evidence "
                            "(up to 500 rows); 'expression' is skipped above 2,500 genes."},
            {"name": "genes_offset", "type": "integer", "default": 0, "min": 0,
             "description": "Offset into the member-gene list. Used with include=genes."},
            {"name": "genes_limit", "type": "integer", "default": 100, "min": 1, "max": 500,
             "description": "Member genes per response. Used with include=genes. The largest "
                            "module has 4,691; /link/gene-module pages through membership too."},
        ],
        "returns": "entity",
        "links": ["self", "html", "parent", "genes"],
        "example_url": "/endoderm-perturbseq/api/v1/module/DE-1",
        "errors": [
            {"status": 400, "code": "ambiguous_name",
             "when": "The name exists in more than one collection (e.g. 'unassigned'); pass ?source=."},
            {"status": 404, "code": "module_not_found", "when": "No module with that name."},
        ],
        "notes": [
            "Gene clusters are reported with their public display name (GC1), never the internal "
            "cluster_1 form; the full alias set is in the 'aliases' field.",
        ],
    },
    {
        "id": "submodule",
        "group": "objects",
        "phase": 2,
        "method": "GET",
        "rule": "/submodule/<submodule>",
        "endpoint": "perturbseq_api_v1.submodule",
        "path": "/submodule/{submodule}",
        "summary": "One submodule. Identical in shape to /module, but pinned to the "
                   "hotspot_submodule collection so the name can never be ambiguous.",
        "path_params": [
            {"name": "submodule", "type": "string", "required": True, "example": "DE-1.1",
             "description": "Submodule name, e.g. DE-1.1."},
        ],
        "query_params": [
            {"name": "include", "type": "enum_list", "enum": "module_include", "default": "",
             "description": "Comma-separated optional lists, as for /module."},
            {"name": "genes_offset", "type": "integer", "default": 0, "min": 0,
             "description": "Offset into the member-gene list. Used with include=genes."},
            {"name": "genes_limit", "type": "integer", "default": 100, "min": 1, "max": 500,
             "description": "Member genes per response. Used with include=genes."},
        ],
        "returns": "entity",
        "links": ["self", "html", "parent"],
        "example_url": "/endoderm-perturbseq/api/v1/submodule/DE-1.1",
        "errors": [{"status": 404, "code": "module_not_found", "when": "No such submodule."}],
        "notes": [],
    },
    {
        "id": "go_term",
        "group": "objects",
        "phase": 2,
        "method": "GET",
        "rule": "/go-term/<go_term>",
        "endpoint": "perturbseq_api_v1.go_term",
        "path": "/go-term/{go_term}",
        "summary": "One Gene Ontology term: its definition and the number of genes annotated "
                   "to it. The gene list and enrichment results are opt-in via include=.",
        "path_params": [
            {"name": "go_term", "type": "string", "required": True, "example": "GO:0030183",
             "description": "GO accession. 'GO:0030183', 'GO_0030183' and '0030183' are all accepted."},
        ],
        "query_params": [
            {"name": "include", "type": "enum_list", "enum": "go_include", "default": "",
             "description": "Optional lists: genes (paged by genes_offset/genes_limit), "
                            "module_enrichment, tf_enrichment."},
            {"name": "genes_offset", "type": "integer", "default": 0, "min": 0,
             "description": "Offset into the annotated-gene list. Used with include=genes."},
            {"name": "genes_limit", "type": "integer", "default": 100, "min": 1, "max": 500,
             "description": "Annotated genes per response. Used with include=genes. The largest "
                            "terms have ~12,800."},
        ],
        "returns": "entity",
        "links": ["self", "html"],
        "example_url": "/endoderm-perturbseq/api/v1/go-term/GO:0030183",
        "errors": [
            {"status": 404, "code": "go_term_not_found",
             "when": "No such GO accession. An unparseable accession is treated as missing, "
                     "not as a client error."},
        ],
        "notes": ["n_genes is the exact annotation count, whether or not include=genes is set."],
    },
    {
        "id": "dataset",
        "group": "objects",
        "phase": 2,
        "method": "GET",
        "rule": "/dataset/<namespace>/<int:dataset_id>",
        "endpoint": "perturbseq_api_v1.dataset",
        "path": "/dataset/{namespace}/{dataset_id}",
        "summary": "One source dataset — either a TF binding experiment or a histone-PTM "
                   "experiment. The two are separate namespaces that happen to share the word.",
        "path_params": [
            {"name": "namespace", "type": "enum", "enum": "dataset_namespace", "required": True,
             "example": "tf", "description": "'tf' for TF binding datasets, 'ptm' for histone PTM."},
            {"name": "dataset_id", "type": "integer", "required": True, "example": "1",
             "description": "Numeric dataset id within that namespace."},
        ],
        "query_params": [],
        "returns": "entity",
        "links": ["self", "tf"],
        "example_url": "/endoderm-perturbseq/api/v1/dataset/tf/1",
        "errors": [{"status": 404, "code": "dataset_not_found", "when": "No dataset with that id."}],
        "notes": [
            "Peak lists are not served. A TF dataset averages ~6,900 peaks over a 121M-row table; "
            "use the original accession in the 'dataset' field to fetch peaks from ENCODE or ReMap.",
        ],
    },
    {
        "id": "atac_peak",
        "group": "objects",
        "phase": 2,
        "method": "GET",
        "rule": "/atac-peak/<peak>",
        "endpoint": "perturbseq_api_v1.atac_peak",
        "path": "/atac-peak/{peak}",
        "summary": "One ATAC peak: coordinates, accessibility across the time course, linked "
                   "genes, and overlapping histone PTM datasets.",
        "path_params": [
            {"name": "peak", "type": "string", "required": True, "example": "58938",
             "description": "Numeric atac_peak_id or a peak name such as ESC_DE_peak_99081."},
        ],
        "query_params": [],
        "returns": "entity",
        "links": ["self", "tfs"],
        "example_url": "/endoderm-perturbseq/api/v1/atac-peak/58938",
        "errors": [{"status": 404, "code": "peak_not_found", "when": "No peak with that id or name."}],
        "notes": [
            "Only the count of TF overlaps is returned here — a peak can overlap 10,000+ TF peaks. "
            "Use /atac-peak/{peak}/tfs for the list.",
        ],
    },
    {
        "id": "atac_peak_tfs",
        "group": "objects",
        "phase": 2,
        "method": "GET",
        "rule": "/atac-peak/<peak>/tfs",
        "endpoint": "perturbseq_api_v1.atac_peak_tfs",
        "path": "/atac-peak/{peak}/tfs",
        "summary": "Transcription factors with a binding peak overlapping this ATAC peak, "
                   "aggregated per TF.",
        "path_params": [
            {"name": "peak", "type": "string", "required": True, "example": "58938",
             "description": "Numeric atac_peak_id or a peak name."},
        ],
        "query_params": [P_PAGE, P_PER_PAGE],
        "returns": "collection",
        "example_url": "/endoderm-perturbseq/api/v1/atac-peak/58938/tfs",
        "errors": [{"status": 404, "code": "peak_not_found", "when": "No peak with that id or name."}],
        "notes": ["total is null: counting distinct TFs on a busy peak is as expensive as the "
                  "query itself, so it is deliberately not computed."],
    },
    {
        "id": "link_tf_gene",
        "group": "links",
        "phase": 2,
        "method": "GET",
        "rule": "/link/tf-gene/<tf>/<gene>",
        "endpoint": "perturbseq_api_v1.link_tf_gene",
        "path": "/link/tf-gene/{tf}/{gene}",
        "summary": "Evidence that a TF regulates a gene: how many datasets place the TF at the "
                   "gene and the strongest binding score, plus the perturbation effect of "
                   "knocking the TF out. Per-dataset and per-peak detail are opt-in.",
        "path_params": [
            {"name": "tf", "type": "string", "required": True, "example": "FOXA2",
             "description": "TF gene symbol."},
            {"name": "gene", "type": "string", "required": True, "example": "SOX17",
             "description": "Target gene symbol."},
        ],
        "query_params": [
            {"name": "include", "type": "enum_list", "enum": "tf_gene_include", "default": "",
             "description": "'datasets' lists each binding dataset with its scores (up to 500; "
                            "CTCF has 1,000+). 'elements' adds peak-level evidence — the most "
                            "expensive query in the API, skipped entirely when there is no "
                            "dataset-level binding to explain."},
        ],
        "returns": "entity",
        "links": ["self", "tf", "gene", "html"],
        "example_url": "/endoderm-perturbseq/api/v1/link/tf-gene/FOXA2/SOX17",
        "errors": [{"status": 404, "code": "gene_not_found", "when": "Either name does not resolve."}],
        "notes": ["binding_evidence is false when no dataset places this TF at this gene; the "
                  "elements block is then empty regardless of include."],
    },
    {
        "id": "link_tf_module",
        "group": "links",
        "phase": 2,
        "method": "GET",
        "rule": "/link/tf-module/<tf>/<module>",
        "endpoint": "perturbseq_api_v1.link_tf_module",
        "path": "/link/tf-module/{tf}/{module}",
        "summary": "Evidence that a TF regulates a module, from both the perturbation screen "
                   "and binding enrichment, with the per-gRNA detail behind it.",
        "path_params": [
            {"name": "tf", "type": "string", "required": True, "example": "ARID1A",
             "description": "TF gene symbol."},
            {"name": "module", "type": "string", "required": True, "example": "DE-1",
             "description": "Module name; the same aliases as /module are accepted."},
        ],
        "query_params": [P_MODULE_SOURCE],
        "returns": "entity",
        "links": ["self", "tf", "module"],
        "example_url": "/endoderm-perturbseq/api/v1/link/tf-module/ARID1A/DE-1",
        "errors": [
            {"status": 400, "code": "ambiguous_name", "when": "The module name needs ?source=."},
            {"status": 404, "code": "module_not_found", "when": "No such module."},
            {"status": 404, "code": "gene_not_found", "when": "No such TF."},
        ],
        "notes": [],
    },
    {
        "id": "link_gene_module",
        "group": "links",
        "phase": 2,
        "method": "GET",
        "rule": "/link/gene-module",
        "endpoint": "perturbseq_api_v1.link_gene_module",
        "path": "/link/gene-module",
        "summary": "Gene-to-module membership edges. At least one filter is required.",
        "path_params": [],
        "query_params": [
            {"name": "gene", "type": "string", "default": "",
             "description": "Return this gene's memberships."},
            {"name": "module", "type": "string", "default": "",
             "description": "Return this module's members."},
            P_MODULE_SOURCE,
            {"name": "include_unassigned", "type": "boolean", "default": False,
             "description": "Include the 'unassigned' pseudo-modules, which hold 5,805 genes "
                            "between them and are a non-result rather than a module."},
            P_PAGE, P_PER_PAGE,
        ],
        "returns": "collection",
        "example_url": "/endoderm-perturbseq/api/v1/link/gene-module?module=DE-1",
        "errors": [
            {"status": 400, "code": "filter_required",
             "when": "None of gene, module or source was supplied — an unfiltered dump is 26K rows."},
        ],
        "notes": [],
    },
    {
        "id": "edges",
        "group": "links",
        "phase": 2,
        "method": "GET",
        "rule": "/edges",
        "endpoint": "perturbseq_api_v1.edges",
        "path": "/edges",
        "summary": "The TF-to-module regulatory network as an edge list, merging perturbation "
                   "and binding evidence.",
        "path_params": [],
        "query_params": [
            {"name": "level", "type": "enum", "enum": "module_source",
             "default": "hotspot_supermodule",
             "description": "Which module collection to build the network over."},
            {"name": "evidence", "type": "enum", "enum": "edge_evidence", "default": "any",
             "description": "Restrict to edges with this kind of support."},
            {"name": "min_odds_ratio", "type": "number", "default": 1.0, "min": 0,
             "description": "Minimum binding enrichment odds ratio."},
            {"name": "max_padj", "type": "number", "default": 0.05, "min": 0, "max": 1,
             "description": "Maximum adjusted p-value for the binding arm."},
            {"name": "min_abs_nes", "type": "number", "default": 0, "min": 0,
             "description": "Minimum |NES| for the perturbation arm."},
            P_PAGE, P_PER_PAGE,
        ],
        "returns": "collection",
        "example_url": "/endoderm-perturbseq/api/v1/edges",
        "errors": [],
        "notes": [
            "At the default supermodule level the whole network is small enough to page through "
            "in full. The submodule level is ~20x larger, so the significance defaults matter more.",
        ],
    },
    {
        "id": "genes",
        "group": "collections",
        "phase": 2,
        "method": "GET",
        "rule": "/genes",
        "endpoint": "perturbseq_api_v1.genes",
        "path": "/genes",
        "summary": "All genes in the dataset, filterable by biotype, chromosome and "
                   "perturbation-library membership.",
        "path_params": [],
        "query_params": [
            {"name": "biotype", "type": "string", "default": "",
             "description": "Exact gene_biotype, e.g. protein_coding or lncRNA."},
            {"name": "chr", "type": "string", "default": "",
             "description": "Chromosome. '7' and 'chr7' are both accepted."},
            {"name": "perturbed", "type": "boolean", "default": None,
             "description": "Restrict to genes that are (or are not) in the perturbation library."},
            P_PAGE, P_PER_PAGE,
        ],
        "returns": "collection",
        "example_url": "/endoderm-perturbseq/api/v1/genes?per_page=5",
        "errors": [],
        "notes": ["Only primary_gene=1 rows are listed, so each symbol appears once."],
    },
    {
        "id": "tfs",
        "group": "collections",
        "phase": 2,
        "method": "GET",
        "rule": "/tfs",
        "endpoint": "perturbseq_api_v1.tfs",
        "path": "/tfs",
        "summary": "Transcription factors. There is no TF table in the database, so which "
                   "population you get is an explicit choice.",
        "path_params": [],
        "query_params": [
            {"name": "set", "type": "enum", "enum": "tf_set", "default": "perturbed",
             "description": "'perturbed' = the 75 TFs targeted in the CRISPR screen; "
                            "'binding' = the ~1,705 genes with a binding dataset; 'all' = either."},
            P_PAGE, P_PER_PAGE,
        ],
        "returns": "collection",
        "example_url": "/endoderm-perturbseq/api/v1/tfs",
        "errors": [],
        "notes": ["A few perturbed TFs carry stale HGNC symbols with no gene_table row; they are "
                  "returned with gene_id null and symbol_status 'stale' rather than dropped."],
    },
    {
        "id": "modules",
        "group": "collections",
        "phase": 2,
        "method": "GET",
        "rule": "/modules",
        "endpoint": "perturbseq_api_v1.modules",
        "path": "/modules",
        "summary": "All modules across the three collections, with sizes and titles. Full "
                   "descriptions are on /module/{module}.",
        "path_params": [],
        "query_params": [
            P_MODULE_SOURCE,
            {"name": "include_unassigned", "type": "boolean", "default": False,
             "description": "Include the 'unassigned' pseudo-modules."},
            P_PAGE, P_PER_PAGE,
        ],
        "returns": "collection",
        "example_url": "/endoderm-perturbseq/api/v1/modules",
        "errors": [],
        "notes": ["Gene clusters carry both display_name (GC1) and the full alias list."],
    },
    {
        "id": "search",
        "group": "collections",
        "phase": 2,
        "method": "GET",
        "rule": "/search",
        "endpoint": "perturbseq_api_v1.search",
        "path": "/search",
        "summary": "Prefix search across genes, synonyms, TFs, modules and GO terms.",
        "path_params": [],
        "query_params": [
            {"name": "q", "type": "string", "required": True, "default": "",
             "description": "Search string, 2-64 characters. Matched as a prefix for names and "
                            "as a substring for GO term descriptions. Wildcards are escaped, "
                            "not interpreted."},
            {"name": "type", "type": "enum_list", "enum": "search_type", "default": "",
             "description": "Restrict to these result types."},
            {"name": "limit", "type": "integer", "default": 25, "min": 1, "max": 100,
             "description": "Maximum results across all types."},
        ],
        "returns": "collection",
        "example_url": "/endoderm-perturbseq/api/v1/search?q=SOX",
        "errors": [
            {"status": 400, "code": "invalid_param",
             "when": "q is missing, shorter than 2 characters, or longer than 64."},
        ],
        "notes": ["Each result carries matched_on, so an exact symbol hit can be told apart from "
                  "a synonym or a GO description match."],
    },
]
