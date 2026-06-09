"""One-time script to add missing indexes to tf-perturbseq-v5.db.

Run once before starting the web server:
    python add_indexes.py

Safe to re-run — all statements use CREATE INDEX IF NOT EXISTS.
Index creation on large tables (100M+ rows) takes several minutes.
"""
import sqlite3
import pathlib
import time

DB = str(pathlib.Path(__file__).parent / "data" / "db" / "tf-perturbseq-v5.db")

INDEXES = [
    # Enables WHERE tf_gene_name=? filter in _query_tf_gene_link and _query_tf_binding_peaks
    # without a full scan of tf_dataset_table (17k rows, but hit on every TF page load)
    ("idx_tf_dataset_tf_gene_name",
     "CREATE INDEX IF NOT EXISTS idx_tf_dataset_tf_gene_name ON tf_dataset_table(tf_gene_name)"),

    # Composite index for JOIN condition: atac_peak_id AND gene_id
    # Used in _query_tf_gene_link UNION subquery and EXISTS checks
    ("idx_atac_tss_peak_gene",
     "CREATE INDEX IF NOT EXISTS idx_atac_tss_peak_gene ON atac_tss_links(atac_peak_id, gene_id)"),

    # Same for multiome_atac_overlaps
    ("idx_multiome_peak_gene",
     "CREATE INDEX IF NOT EXISTS idx_multiome_peak_gene ON multiome_atac_overlaps(atac_peak_id, gene_id)"),
]

db = sqlite3.connect(DB)
for name, sql in INDEXES:
    print(f"Creating {name} ...", end=" ", flush=True)
    t0 = time.time()
    db.execute(sql)
    db.commit()
    print(f"done ({time.time() - t0:.1f}s)")

db.close()
print("All indexes created.")
