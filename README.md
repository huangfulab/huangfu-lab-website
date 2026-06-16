# Hotspot Network Viewer — Webapp

Interactive browser for Perturb-seq co-expression modules derived from [Hotspot](https://github.com/yoseflab/Hotspot) local-correlation z-scores in an ESC→DE (embryonic stem cell to definitive endoderm) differentiation dataset.

## Quick start

```bash
pip install -r ../../requirements.txt
python app.py
# Open http://localhost:5050
```

> Must be run from this directory (`pipeline/webapp/`) because blueprints resolve template/static paths relative to `__file__`. Alternatively: `python pipeline/webapp/app.py` from the repo root.

---

## App structure

`app.py` is the entry point. It registers two blueprints and serves the network view directly:

| Mount | File | Description |
|---|---|---|
| `/tf-perturbseq/` | `perturbseq_bp.py` | Main biology browser: TF, module/submodule, gene, GO term, regulatory-element, and TF–gene link pages |
| `/timecourse/` | `timecourse_bp.py` | ESC→DE time-course network (gene clusters × peak clusters) |
| `/network` | `app.py` | Cytoscape interactive network of supermodule/submodule co-expression |

### Databases

Both blueprints need SQLite databases built before first run:

| Database | Location | Built by |
|---|---|---|
| `perturbseq.db` | `../../data/perturbseq.db` | external pipeline |
| `timecourse.db` | `timecourse.db` (this dir) | `python ../../pipeline/preprocess_timecourse.py` |

`perturbseq_bp.py` opens both databases — gene-cluster (GC) and peak-cluster (PC) pages live under the `/tf-perturbseq/` URL scheme (controlled by `PERTURBSEQ_PREFIX` in `perturbseq_bp.py`) but pull timecourse data.

---

## Data hierarchy

**Co-expression modules** (two-tier, from Hotspot z-scores):

- **Supermodules** — 15 broad programs (`DE-0` … `DE-14`)
- **Submodules** — ~200 finer clusters (e.g. `DE-6.5`); name encodes parent supermodule

**Timecourse clusters** (ESC→DE time series):

- **Gene clusters** — GC1–GC6 (mfuzz temporal expression clusters)
- **Peak clusters** — PC1–PC4, PC6–PC7 (ATAC-seq temporal clusters)

**TF regulatory evidence** (two orthogonal sources):

- **Perturbation** (Perturb-seq GSEA) — TF knockdown → module activity shift (NES, padj)
- **Binding** (ChIP-seq / TOBIAS footprints) — TF peak enrichment in ATAC clusters (odds ratio, FDR)

---

## Templates and static assets

```
templates/
  landing.html          # / landing page
  index.html            # /network Cytoscape view
  modules.html          # /modules overview
  perturbseq/           # /tf-perturbseq/* pages (tf, module, submodule, gene, go_term, element, ...)
  timecourse/           # /timecourse/* pages (index, cluster, dual, sankey, ...)

static/
  js/                   # Custom JS (tf-network.js, etc.)
  css/                  # Stylesheets
  cytoscape*.js         # Cytoscape + layout plugins (bundled)
  vis-network.min.js    # Vis.js network (bundled)
```

---

## Pipeline scripts (run from repo root)

| Script | Purpose |
|---|---|
| `pipeline/preprocess_timecourse.py` | Build `timecourse.db` from `data/04-timecourse_data/` |
| `pipeline/generate_submodule_descriptions.py` | Generate LLM functional descriptions for submodules |
| `pipeline/collate_submodule_summaries.py` | Collate submodule gene lists into JSON for LLM input |
| `pipeline/extract_mfuzz_profiles.R` | Extract mfuzz soft-clustering profiles to TSV |

---

## Key source data files

| Path | Description |
|---|---|
| `data/01-hotspot/DE_hotspot_pca100_gemgroup-local_correlation_z.parquet` | 5,895 × 5,895 gene–gene co-expression z-score matrix |
| `data/tables/` | Source TSVs for most `perturbseq.db` tables |
| `data/datasets/` | Module assignments, pathway enrichment, TF enrichment TSVs |
| `data/07-2026_05_26_data/` | TF–gene genomic mapping (ChIP-seq peaks, ATAC peaks, junction tables) |
| `data/06-llm_submodule_descriptions/submodule_descriptions-v3.tsv` | LLM-generated submodule descriptions |
| `gene_name_map.tsv` | Gene symbol → Ensembl ID map (links Perturb-seq genes to timecourse cluster members) |
| `networks/` | Node/edge TSVs for the `/network` Cytoscape view |
