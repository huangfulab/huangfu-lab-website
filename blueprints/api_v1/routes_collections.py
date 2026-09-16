"""Collection and search endpoints.

Every listing orders by a unique column. Ordering by a non-unique one would let
rows shift between pages on a database this size, silently duplicating and
dropping records for anyone paging through.
"""
from .catalog import ENUMS
from .core import (
    QUERY_DEADLINE_COLLECTION,
    ApiError,
    api_db,
    api_v1_bp,
    arg_bool,
    arg_enum,
    arg_enum_list,
    arg_int,
    arg_str,
    collection,
    html_link,
    page_params,
    path_link,
    query,
    reject_unknown_args,
    scalar,
)
from .names import (
    HIDDEN_MODULES,
    like_escape,
    module_aliases,
    module_display_name,
    norm_chr,
    prefix_range,
)

SEARCH_MIN_Q = 2
SEARCH_MAX_Q = 64

# NULL-guarded filters rather than a WHERE clause assembled at request time, so
# no user value can reach the SQL text. None of these columns is indexed, so
# there is no selectivity to give up.
_GENE_FILTERS = (
    "  AND (? IS NULL OR g.gene_biotype = ?) "
    "  AND (? IS NULL OR g.chr = ?) "
    "  AND (? IS NULL OR g.in_perturbation_library = ?) "
)


@api_v1_bp.route("/genes")
def genes():
    reject_unknown_args("biotype", "chr", "perturbed", "page", "per_page")
    biotype = arg_str("biotype")
    chrom = norm_chr(arg_str("chr"))
    perturbed = arg_bool("perturbed")
    page, per_page, offset = page_params()

    db = api_db(QUERY_DEADLINE_COLLECTION)
    biotype = biotype or None
    chrom = chrom or None
    perturbed = None if perturbed is None else (1 if perturbed else 0)
    params = (biotype, biotype, chrom, chrom, perturbed, perturbed)

    total = scalar(
        db,
        "SELECT COUNT(*) FROM gene_table g WHERE g.primary_gene = 1 " + _GENE_FILTERS,
        params,
        default=0,
    )
    rows = query(
        db,
        "SELECT g.gene_id, g.ensg_id, g.gene_name, g.gene_biotype, g.chr, "
        "       g.chrom_start, g.chrom_end, g.description, g.in_perturbation_library, "
        "       pc.publication_count "
        "FROM gene_table g "
        "LEFT JOIN gene_publication_count pc ON pc.gene_id = g.gene_id "
        "WHERE g.primary_gene = 1 " + _GENE_FILTERS +
        "ORDER BY g.gene_name, g.gene_id LIMIT ? OFFSET ?",
        (*params, per_page, offset),
    )
    for r in rows:
        r["id"] = r["gene_name"]
        r["type"] = "gene"
        r["in_perturbation_library"] = bool(r["in_perturbation_library"])
        r["url"] = path_link("gene", r["gene_name"])
    return collection(rows, page=page, per_page=per_page, total=total)


@api_v1_bp.route("/tfs")
def tfs():
    reject_unknown_args("set", "page", "per_page")
    tf_set = arg_enum("set", ENUMS["tf_set"], default="perturbed")
    page, per_page, offset = page_params()

    db = api_db(QUERY_DEADLINE_COLLECTION)

    # Perturbed and binding TFs are different populations drawn from different
    # tables; 'all' is their union, keyed on symbol because a few perturbed
    # entries carry stale HGNC symbols with no gene_table row.
    perturbed = {
        r["gene_name"]: r
        for r in query(
            db,
            "SELECT s.gene_name, s.gene_id, COUNT(DISTINCT s.module) AS n_modules_regulated "
            "FROM gsea_tf_table s WHERE s.gene_set_collection = 'DE-hotspot_modules' "
            "GROUP BY s.gene_name, s.gene_id",
        )
    } if tf_set in ("perturbed", "all") else {}

    binding = {
        r["gene_name"]: r
        for r in query(
            db,
            "SELECT g.gene_id, g.gene_name, COUNT(DISTINCT tdg.dataset_id) AS n_binding_datasets "
            "FROM tf_dataset_gene_table tdg JOIN gene_table g ON g.gene_id = tdg.gene_id "
            "GROUP BY g.gene_id",
        )
    } if tf_set in ("binding", "all") else {}

    merged = {}
    for name, r in perturbed.items():
        merged[name] = {
            "id": name, "type": "tf", "gene_name": name, "gene_id": r["gene_id"],
            "is_perturbed_tf": True, "is_binding_tf": False,
            "n_modules_regulated": r["n_modules_regulated"], "n_binding_datasets": 0,
            "symbol_status": "stale" if r["gene_id"] is None else "current",
        }
    for name, r in binding.items():
        entry = merged.setdefault(name, {
            "id": name, "type": "tf", "gene_name": name, "gene_id": r["gene_id"],
            "is_perturbed_tf": False, "is_binding_tf": False,
            "n_modules_regulated": 0, "n_binding_datasets": 0, "symbol_status": "current",
        })
        entry["is_binding_tf"] = True
        entry["n_binding_datasets"] = r["n_binding_datasets"]
        if entry["gene_id"] is None:
            entry["gene_id"] = r["gene_id"]

    ordered = sorted(merged.values(), key=lambda r: r["gene_name"])
    for r in ordered:
        r["url"] = path_link("tf", r["gene_name"])
    return collection(
        ordered[offset:offset + per_page], page=page, per_page=per_page, total=len(ordered)
    )


@api_v1_bp.route("/modules")
def modules():
    reject_unknown_args("source", "include_unassigned", "page", "per_page")
    source = arg_enum("source", ENUMS["module_source"])
    include_unassigned = arg_bool("include_unassigned", False)
    page, per_page, offset = page_params()

    db = api_db(QUERY_DEADLINE_COLLECTION)
    rows = query(
        db,
        # Title only: the paragraph-length descriptions are served by /module/{name}
        # and would triple the size of every listing page.
        "SELECT m.module_id, m.module_name, m.source, m.size, d.title, "
        "       (SELECT COUNT(*) FROM gene_module_table gm WHERE gm.module_id = m.module_id) "
        "         AS n_genes "
        "FROM module_table m "
        "LEFT JOIN module_description d ON d.module_id = m.module_id "
        "WHERE (? IS NULL OR m.source = ?) "
        "  AND (? = 1 OR m.module_name <> 'unassigned') "
        "ORDER BY m.source, m.module_name, m.module_id",
        (source, source, 1 if include_unassigned else 0),
    )
    visible = [r for r in rows if (r["source"], r["module_name"]) not in HIDDEN_MODULES]
    for r in visible:
        r["display_name"] = module_display_name(r["module_name"], r["source"])
        r["aliases"] = module_aliases(r["module_name"], r["source"])
        r["id"] = r["display_name"]
        r["type"] = "module"
        r["url"] = path_link("module", r["display_name"])
    return collection(
        visible[offset:offset + per_page], page=page, per_page=per_page, total=len(visible)
    )


def _search_genes(db, q, limit):
    """Prefix match via a range predicate rather than LIKE.

    Gene symbols are uppercase in all but 204 of 24,960 rows, so probing both the
    raw and upper-cased forms covers the table without falling back to a scan.
    """
    out = {}
    for variant in {q, q.upper()}:
        bounds = prefix_range(variant)
        if not bounds:
            continue
        for r in query(
            db,
            "SELECT gene_id, gene_name, gene_biotype, in_perturbation_library "
            "FROM gene_table WHERE gene_name >= ? AND gene_name < ? AND primary_gene = 1 "
            "ORDER BY gene_name LIMIT ?",
            (*bounds, limit),
        ):
            out.setdefault(r["gene_name"], r)
    return list(out.values())


def _search_synonyms(db, q, limit):
    out = {}
    for variant in {q, q.upper()}:
        bounds = prefix_range(variant)
        if not bounds:
            continue
        for r in query(
            db,
            "SELECT s.synonym, g.gene_id, g.gene_name FROM gene_synonym s "
            "JOIN gene_table g ON g.gene_id = s.gene_id "
            "WHERE s.synonym >= ? AND s.synonym < ? ORDER BY s.synonym LIMIT ?",
            (*bounds, limit),
        ):
            out.setdefault(r["synonym"], r)
    return list(out.values())


@api_v1_bp.route("/search")
def search():
    reject_unknown_args("q", "type", "limit")
    q = arg_str("q", max_len=SEARCH_MAX_Q + 1)
    types = arg_enum_list("type", ENUMS["search_type"]) or set(ENUMS["search_type"])
    limit = arg_int("limit", 25, minimum=1, maximum=100)

    if len(q) < SEARCH_MIN_Q:
        raise ApiError(400, f"Parameter 'q' must be at least {SEARCH_MIN_Q} characters.")
    if len(q) > SEARCH_MAX_Q:
        raise ApiError(400, f"Parameter 'q' must be at most {SEARCH_MAX_Q} characters.")

    db = api_db(QUERY_DEADLINE_COLLECTION)
    results = []

    if {"gene", "tf"} & types:
        for r in _search_genes(db, q, limit):
            is_tf = bool(r["in_perturbation_library"])
            kind = "tf" if is_tf else "gene"
            if kind in types:
                results.append({
                    "type": kind, "id": r["gene_name"], "name": r["gene_name"],
                    "display_name": r["gene_name"], "gene_id": r["gene_id"],
                    "matched_on": "gene_name", "url": path_link(kind, r["gene_name"]),
                    "html_url": html_link("gene", r["gene_name"]),
                })

    if "synonym" in types:
        known = {r["id"] for r in results}
        for r in _search_synonyms(db, q, limit):
            if r["gene_name"] in known:
                continue
            results.append({
                "type": "synonym", "id": r["gene_name"], "name": r["synonym"],
                "display_name": f"{r['synonym']} → {r['gene_name']}", "gene_id": r["gene_id"],
                "matched_on": "synonym", "url": path_link("gene", r["gene_name"]),
                "html_url": html_link("gene", r["gene_name"]),
            })

    if {"module", "submodule", "gene_cluster"} & types:
        # 313 rows, so an escaped LIKE is cheap here; the escape still matters so
        # a literal % or _ in the query is not treated as a wildcard.
        for r in query(
            db,
            "SELECT module_id, module_name, source, size FROM module_table "
            "WHERE module_name LIKE ? ESCAPE '\\' ORDER BY source, module_name LIMIT ?",
            (like_escape(q) + "%", limit),
        ):
            if (r["source"], r["module_name"]) in HIDDEN_MODULES:
                continue
            kind = {"hotspot_supermodule": "module", "hotspot_submodule": "submodule",
                    "mfuzz_k7": "gene_cluster"}[r["source"]]
            if kind not in types:
                continue
            display = module_display_name(r["module_name"], r["source"])
            results.append({
                "type": kind, "id": display, "name": display, "display_name": display,
                "module_id": r["module_id"], "source": r["source"], "matched_on": "module_name",
                "url": path_link("module", display), "html_url": html_link("module", display),
            })

    if "gene_cluster" in types and q.upper().startswith(("GC", "TC-")):
        digits = q.upper().removeprefix("GC").removeprefix("TC-")
        if digits.isdigit() and not any(r["id"] == f"GC{int(digits)}" for r in results):
            row = query(
                db,
                "SELECT module_id, module_name, source FROM module_table "
                "WHERE module_name = ? AND source = 'mfuzz_k7'",
                (f"cluster_{int(digits)}",),
            )
            for r in row:
                if (r["source"], r["module_name"]) in HIDDEN_MODULES:
                    continue
                display = module_display_name(r["module_name"], r["source"])
                results.append({
                    "type": "gene_cluster", "id": display, "name": display,
                    "display_name": display, "module_id": r["module_id"], "source": r["source"],
                    "matched_on": "alias", "url": path_link("module", display),
                    "html_url": html_link("module", display),
                })

    if "go_term" in types:
        if q.upper().startswith("GO:"):
            go_rows = query(
                db,
                "SELECT go_id, go_accession, go_name, namespace FROM go_term_table "
                "WHERE go_accession = ? AND is_obsolete = 0",
                (q.upper(),),
            )
            matched = "go_accession"
        else:
            # Substring search over 48K rows is an unavoidable scan, but it is the
            # only scan in the request and measures ~45 ms.
            go_rows = query(
                db,
                "SELECT go_id, go_accession, go_name, namespace FROM go_term_table "
                "WHERE is_obsolete = 0 AND go_name LIKE ? ESCAPE '\\' "
                "ORDER BY LENGTH(go_name) LIMIT ?",
                ("%" + like_escape(q) + "%", limit),
            )
            matched = "go_name"
        for r in go_rows:
            results.append({
                "type": "go_term", "id": r["go_accession"], "name": r["go_name"],
                "display_name": f"{r['go_name']} ({r['go_accession']})",
                "namespace": r["namespace"], "matched_on": matched,
                "url": path_link("go-term", r["go_accession"]),
                # The site's canonical GO page is "<name> (<accession>)"; the
                # accession alone answers with a 301.
                "html_url": html_link("go", f"{r['go_name']} ({r['go_accession']})"),
            })

    exact = q.upper()
    results.sort(key=lambda r: (r["name"].upper() != exact, len(r["name"]), r["name"]))
    trimmed = results[:limit]
    return collection(trimmed, page=1, per_page=limit, total=len(trimmed), paginated=False)
