"""Entity name resolution — the single source of aliasing truth for the API.

Every rule about how a user-supplied name maps onto the database lives here, so
there is one place to check when a lookup returns something surprising.
"""
import re

from .core import ApiError, clean_name, query, query_one

# The DB stores mfuzz gene clusters as cluster_N. The public display name is GCN,
# the legacy internal id is gene_cluster_N, and the legacy URL form was TC-N.
# "GC" never appears in the database.
_MFUZZ_RE = re.compile(r"^(?:GC|TC-|cluster_|gene_cluster_)(\d+)$", re.IGNORECASE)
_SUBMODULE_RE = re.compile(r"^DE-\d+\.\d+$")
_SUPERMODULE_RE = re.compile(r"^DE-\d+$")

MFUZZ_SOURCE = "mfuzz_k7"

# Hidden in the browser UI (HIDDEN_TC in perturbseq_bp.py); hidden here too, so
# the API and the site agree. peak_cluster_5 has no module_table row.
HIDDEN_MODULES = {(MFUZZ_SOURCE, "cluster_7")}

# Note there is deliberately no gene_set_collection map for tf_module_enrichment:
# that table labels the mfuzz collection 'mfuzz' while module_table calls it
# 'mfuzz_k7', so every query here drives off module_id instead and never has to
# cross the two vocabularies. Filtering it by gene_set_collection would also turn
# an indexed seek into a full scan of 521K rows.

# gsea_tf_table holds 28 distinct gene_set_collection values, 20 of them
# near-identical mfuzz variants. The app canonicalises on these two.
GSEA_COLLECTION = {
    "hotspot_supermodule": "DE-hotspot_modules",
    "hotspot_submodule": "DE-hotspot_modules",
    MFUZZ_SOURCE: "ESC_DE-gene_clustering_data_var0.3_k7_top1000_DE",
}


def normalize_module(name):
    """Return (db_module_name, inferred_source_or_None)."""
    name = clean_name(name, "module")
    mfuzz = _MFUZZ_RE.match(name)
    if mfuzz:
        return f"cluster_{int(mfuzz.group(1))}", MFUZZ_SOURCE
    if _SUBMODULE_RE.match(name):
        return name, "hotspot_submodule"
    if _SUPERMODULE_RE.match(name):
        return name, "hotspot_supermodule"
    return name, None


def module_display_name(module_name, source):
    if source == MFUZZ_SOURCE and module_name.startswith("cluster_"):
        return "GC" + module_name.split("_", 1)[1]
    return module_name


def module_aliases(module_name, source):
    if source != MFUZZ_SOURCE or not module_name.startswith("cluster_"):
        return []
    n = module_name.split("_", 1)[1]
    return [f"GC{n}", f"TC-{n}", f"gene_cluster_{n}", module_name]


def resolve_module(db, name, source=None, include_hidden=False):
    """Resolve a module name to its module_table row.

    Raises 400 when the name is ambiguous across collections ('unassigned' exists
    as both a supermodule and a submodule) and 404 when it does not exist.
    """
    db_name, inferred = normalize_module(name)
    source = source or inferred

    if source:
        rows = query(
            db,
            "SELECT module_id, module_name, source, size FROM module_table "
            "WHERE module_name = ? AND source = ?",
            (db_name, source),
        )
    else:
        rows = query(
            db,
            "SELECT module_id, module_name, source, size FROM module_table "
            "WHERE module_name = ?",
            (db_name,),
        )

    if not rows:
        raise ApiError(404, f"No module named {name!r}.", "module_not_found")
    if len(rows) > 1:
        sources = ", ".join(sorted(r["source"] for r in rows))
        raise ApiError(
            400,
            f"Module name {name!r} exists in more than one collection ({sources}). "
            f"Add ?source= to disambiguate.",
            "ambiguous_name",
        )

    row = rows[0]
    if not include_hidden and (row["source"], row["module_name"]) in HIDDEN_MODULES:
        raise ApiError(404, f"No module named {name!r}.", "module_not_found")

    row["display_name"] = module_display_name(row["module_name"], row["source"])
    row["aliases"] = module_aliases(row["module_name"], row["source"])
    return row


def resolve_gene(db, name):
    """Resolve a symbol, Ensembl ID or synonym to a gene_table row.

    Prefers primary_gene=1 — 129 rows share a symbol with a primary entry.
    Reports what was matched so the caller can surface it in meta.
    """
    name = clean_name(name, "gene")

    row = query_one(
        db,
        "SELECT gene_id, gene_name FROM gene_table WHERE gene_name = ? AND primary_gene = 1 LIMIT 1",
        (name,),
    )
    matched_on = "symbol"

    if not row:
        row = query_one(
            db, "SELECT gene_id, gene_name FROM gene_table WHERE gene_name = ? LIMIT 1", (name,)
        )
    if not row:
        row = query_one(
            db, "SELECT gene_id, gene_name FROM gene_table WHERE ensg_id = ? LIMIT 1", (name,)
        )
        matched_on = "ensembl_id" if row else matched_on
    if not row:
        row = query_one(
            db,
            "SELECT g.gene_id, g.gene_name FROM gene_synonym s "
            "JOIN gene_table g ON g.gene_id = s.gene_id WHERE s.synonym = ? "
            "ORDER BY g.primary_gene DESC LIMIT 1",
            (name,),
        )
        matched_on = "synonym" if row else matched_on

    if not row:
        raise ApiError(
            404,
            f"No gene matching {name!r}. Symbols, Ensembl gene IDs and synonyms are accepted.",
            "gene_not_found",
        )

    row["matched_on"] = matched_on
    row["resolved_from"] = name if name != row["gene_name"] else None
    return row


def gene_meta(gene_row):
    meta = {}
    if gene_row.get("resolved_from"):
        meta["resolved_from"] = gene_row["resolved_from"]
        meta["match_type"] = gene_row["matched_on"]
    return meta


def norm_chr(value):
    """Chromosomes are stored bare ('1', 'X'), but users will type 'chr7'."""
    value = (value or "").strip()
    return value[3:] if value.lower().startswith("chr") else value


def like_escape(value):
    r"""Escape LIKE metacharacters. Without this, q='100%' silently becomes a
    wildcard search and q='%' matches everything."""
    return value.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def prefix_range(value):
    """Half-open [lo, hi) bounds for a prefix match.

    SQLite only uses an index for LIKE when case_sensitive_like is ON or the
    index is NOCASE — neither holds here, so every LIKE 'X%' is a full scan.
    A range predicate on the same column is index-served instead.
    """
    if not value:
        return None
    return value, value[:-1] + chr(ord(value[-1]) + 1)


def resolve_go_accession(value):
    """Accept 'GO:0030183', 'GO_0030183' or '0030183'.

    Anything unrecognised is passed through untouched so the lookup simply finds
    nothing and the caller returns 404 — a malformed accession is a missing term,
    not a client error worth a different status.
    """
    value = clean_name(value, "GO term").upper().replace("GO_", "GO:")
    if value.isdigit():
        return f"GO:{int(value):07d}"
    return value


def tf_dataset_labels(db, gene_id):
    """Dataset-level TF labels for a gene.

    tf_dataset_table.tf_gene_name holds modified names (GATA3_Nter, POU5F1_M,
    NOTCH1_NICD) that do not equi-join to gene_table, so enrichment lookups keyed
    on that column need the full label set for a gene, not just its symbol.
    """
    return [
        r["tf_gene_name"]
        for r in query(
            db,
            "SELECT DISTINCT td.tf_gene_name FROM tf_dataset_gene_table tdg "
            "JOIN tf_dataset_table td ON td.dataset_id = tdg.dataset_id "
            "WHERE tdg.gene_id = ? ORDER BY td.tf_gene_name",
            (gene_id,),
        )
    ]
