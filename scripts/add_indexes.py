"""Add missing indexes to the active Perturb-Seq database.

    python scripts/add_indexes.py            # cheap tier (seconds to ~2 min)
    python scripts/add_indexes.py --list     # show what would be created, then exit

Safe to re-run: every statement uses CREATE INDEX IF NOT EXISTS.

This writes to a ~33 GB database. Set data/db/db_maintenance.txt to TRUE before
running against a live deployment — the site shows a banner and the public API
returns 503 while it is set.
"""
import argparse
import pathlib
import sqlite3
import sys
import time

DB = str(pathlib.Path(__file__).parent.parent / "data" / "db" / "tf-perturbseq-v6.db")

# Cheap tier: all of these are on tables of 4.1M rows or fewer. Every one backs a
# lookup the public API performs on most requests.
INDEXES = [
    # module_name is the most natural module lookup in the API and was unindexed.
    # Composite with source because the name alone is not unique — 'unassigned'
    # exists as both a supermodule and a submodule.
    ("idx_module_name_source", "module_table",
     "CREATE INDEX IF NOT EXISTS idx_module_name_source ON module_table(module_name, source)"),

    # Backs _query_tf_gene_link / _query_tf_binding_peaks and the API's TF lookups,
    # replacing a 17.6K-row scan hit on every TF page load.
    ("idx_tf_dataset_tf_gene_name", "tf_dataset_table",
     "CREATE INDEX IF NOT EXISTS idx_tf_dataset_tf_gene_name ON tf_dataset_table(tf_gene_name)"),

    # grna_table is indexed on gene_id but not gene_name, which several existing
    # CTEs filter on.
    ("idx_grna_gene_name", "grna_table",
     "CREATE INDEX IF NOT EXISTS idx_grna_gene_name ON grna_table(gene_name)"),

    # The is-this-a-TF test runs on nearly every gene and TF request.
    ("idx_gsea_tf_gene_name_coll", "gsea_tf_table",
     "CREATE INDEX IF NOT EXISTS idx_gsea_tf_gene_name_coll "
     "ON gsea_tf_table(gene_name, gene_set_collection)"),

    # atac_peak_name is unique but unindexed; lets a peak be fetched by name.
    ("idx_atac_peak_name", "atac_peak_table",
     "CREATE INDEX IF NOT EXISTS idx_atac_peak_name ON atac_peak_table(atac_peak_name)"),

    # Makes the module-anchored binding-enrichment query a covering seek.
    ("idx_enrichment_module_tf", "tf_module_enrichment",
     "CREATE INDEX IF NOT EXISTS idx_enrichment_module_tf "
     "ON tf_module_enrichment(module_id, tf_gene_name)"),

    # Composite JOIN keys used by the peak/gene link queries.
    ("idx_atac_tss_peak_gene", "atac_tss_links",
     "CREATE INDEX IF NOT EXISTS idx_atac_tss_peak_gene ON atac_tss_links(atac_peak_id, gene_id)"),
    ("idx_multiome_peak_gene", "multiome_atac_overlaps",
     "CREATE INDEX IF NOT EXISTS idx_multiome_peak_gene "
     "ON multiome_atac_overlaps(atac_peak_id, gene_id)"),

    # 4.1M rows, 1-2 minutes. de_results has separate grna_id and gene_id indexes
    # but no composite, so a single (gRNA, gene) lookup reads every row for the gene.
    ("idx_de_results_gene_grna", "de_results",
     "CREATE INDEX IF NOT EXISTS idx_de_results_gene_grna ON de_results(gene_id, grna_id)"),
]

# Deliberately NOT included: atac_tf_overlaps(atac_peak_id, peak_id, overlap_bp)
# would take 1-3 hours to build and add 3-4 GB, for a ~30% improvement on one
# gated endpoint. Revisit only if that endpoint proves heavily used.


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="show planned indexes and exit")
    args = ap.parse_args()

    if args.list:
        for name, table, _ in INDEXES:
            print(f"{name:32s} on {table}")
        return 0

    if not pathlib.Path(DB).exists():
        print(f"Database not found: {DB}", file=sys.stderr)
        return 2

    db = sqlite3.connect(DB)
    existing = {r[0] for r in db.execute(
        "SELECT name FROM sqlite_master WHERE type='index'")}

    print(f"{DB}\n")
    created = skipped = 0
    for name, table, sql in INDEXES:
        if name in existing:
            print(f"  exists   {name}")
            skipped += 1
            continue
        print(f"  creating {name} on {table} ...", end=" ", flush=True)
        t0 = time.time()
        db.execute(sql)
        db.commit()
        print(f"done ({time.time() - t0:.1f}s)")
        created += 1

    db.close()
    print(f"\n{created} created, {skipped} already present.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
