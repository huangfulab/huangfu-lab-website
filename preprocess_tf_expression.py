#!/usr/bin/env python3
"""
Pre-calculate TFs with max TPM > 1 across any sample and store in perturbseq.db.
Run once (or whenever the expression data changes):
    python pipeline/webapp/preprocess_tf_expression.py
"""
from pathlib import Path
import sqlite3
import csv

ROOT    = Path(__file__).resolve().parent.parent.parent
EXPR    = ROOT / "data" / "03-gene_expression" / "ESC_DE-merged_expression.tsv"
DB_PATH = ROOT / "data" / "perturbseq.db"

# Stream the TSV row-by-row — avoids loading 587 k rows into memory
max_tpm: dict[str, float] = {}
with open(EXPR) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        g, t = row['gene_name'], float(row['tpm'])
        if t > max_tpm.get(g, 0.0):
            max_tpm[g] = t

expressed_genes = {g for g, t in max_tpm.items() if t > 1.0}
print(f"Genes with max TPM > 1 across all samples: {len(expressed_genes)}")

conn = sqlite3.connect(DB_PATH)
tfs_in_db = {r[0] for r in conn.execute(
    "SELECT DISTINCT tf FROM binding_edges "
    "UNION SELECT DISTINCT tf FROM perturbation_edges"
).fetchall()}
print(f"Unique TFs in DB: {len(tfs_in_db)}")

expressed_tfs = expressed_genes & tfs_in_db
print(f"Expressed TFs (overlap with DB): {len(expressed_tfs)}")

conn.execute("DROP TABLE IF EXISTS expressed_tfs")
conn.execute("CREATE TABLE expressed_tfs (tf TEXT PRIMARY KEY)")
conn.executemany("INSERT INTO expressed_tfs VALUES (?)", [(tf,) for tf in sorted(expressed_tfs)])
conn.commit()
conn.close()
print("Done — expressed_tfs table written to perturbseq.db")
