# Migration Query Map: perturbseq.db → TfDatabase

Every SQL query in `perturbseq_bp.py` mapped to its new equivalent under the TfDatabase schema.

**Key constants (TfDatabase):**
- Supermodules: `module_table.source = 'hotspot_supermodule'`
- Submodules: `module_table.source = 'hotspot_submodule'`
- GSEA collection: `gsea_tf_table.gene_set_collection = 'DE-hotspot_modules'`
- GSEA module collection: `gsea_tf_table.module_collection = 'hotspot_supermodule'` or `'hotspot_submodule'`
- Samples: `ES_0h_1..3`, `DE_12h_1..3`, ..., `DE_72h_1..3` — strip `_\d+$` to get timepoint

---

## Module-level caches

### `_load_expressed_tfs()`
```sql
-- OLD (expressed_tfs table removed)
SELECT tf FROM expressed_tfs

-- NEW (TFs with perturbation data in library)
SELECT DISTINCT gene_name FROM gsea_tf_table
-- OR more precisely:
SELECT gene_name FROM gene_table WHERE in_perturbation_library=1
```

### `_load_perturbed_tfs()`
```sql
-- OLD
SELECT DISTINCT tf FROM perturbation_edges

-- NEW
SELECT DISTINCT gene_name FROM gsea_tf_table
```

---

## Navigation lists — `_nav_lists(db)`

```sql
-- OLD: active TFs
SELECT DISTINCT tf FROM tf_descriptions ORDER BY 1

-- NEW
SELECT gene_name FROM gene_table WHERE in_perturbation_library=1 ORDER BY gene_name

-- OLD: modules
SELECT DISTINCT module FROM module_descriptions WHERE module NOT LIKE 'GC%' ORDER BY 1

-- NEW
SELECT module_name FROM module_table
WHERE source='hotspot_supermodule' AND module_name != 'unassigned'
ORDER BY module_name
```

---

## Index page — `index()`

```sql
-- OLD: TF list
SELECT DISTINCT tf FROM perturbation_edges
UNION SELECT DISTINCT tf FROM binding_edges ORDER BY 1

-- NEW
SELECT gene_name FROM gene_table WHERE in_perturbation_library=1 ORDER BY gene_name

-- OLD: module list
SELECT DISTINCT module FROM perturbation_edges WHERE module_type='supermodule'
UNION SELECT DISTINCT module FROM binding_edges WHERE module_type='supermodule' ORDER BY 1

-- NEW
SELECT module_name FROM module_table
WHERE source='hotspot_supermodule' AND module_name != 'unassigned'
ORDER BY module_name
```

---

## Supermodule page — `_super_page(mod_name)`

```sql
-- OLD: existence check
SELECT 1 FROM module_genes WHERE module=? LIMIT 1

-- NEW
SELECT module_id FROM module_table
WHERE module_name=? AND source='hotspot_supermodule' LIMIT 1

-- OLD: description
SELECT * FROM module_descriptions WHERE module=?

-- NEW: table removed — keep submodule_descriptions-v3.tsv externally or add new column

-- OLD: expression (pre-aggregated)
SELECT timepoint, mean_tpm, sd_tpm, n_genes FROM module_expression WHERE module=? ORDER BY ...

-- NEW: call _module_expr_by_timepoint(db, module_id)  [new helper]

-- OLD: member genes
SELECT gene FROM module_genes WHERE module=?

-- NEW
SELECT gene_name FROM gene_module_table WHERE module_id=? AND gene_name IS NOT NULL

-- OLD: submodule listing
SELECT id, n_genes, color, within_mean_z FROM submodule_nodes WHERE supermodule=? ORDER BY n_genes DESC

-- NEW (note: color and within_mean_z not in new schema)
SELECT module_id, module_name, size FROM module_table
WHERE source='hotspot_submodule' AND module_name LIKE ?  -- LIKE = mod_name + '.%'
ORDER BY size DESC

-- OLD: submodule gene lookup
SELECT submodule, gene FROM submodule_genes WHERE submodule IN (...)

-- NEW
SELECT module_id, gene_name FROM gene_module_table
WHERE module_id IN (...) AND gene_name IS NOT NULL
```

---

## Submodule page — `_sub_page(name)`

```sql
-- OLD: node metadata
SELECT * FROM submodule_nodes WHERE id=?

-- NEW (note: color, within_mean_z, supermodule not in schema; derive supermodule from name)
SELECT module_id, module_name, source, size FROM module_table
WHERE module_name=? AND source='hotspot_submodule'

-- OLD: member genes
SELECT gene FROM submodule_genes WHERE submodule=? ORDER BY gene

-- NEW
SELECT gene_name FROM gene_module_table WHERE module_id=? AND gene_name IS NOT NULL ORDER BY gene_name

-- OLD: submodule-submodule edges (REMOVED — table does not exist in TfDatabase)
SELECT source, target, mean_z, n_pairs FROM submodule_edges
WHERE source=? OR target=? ORDER BY ABS(mean_z) DESC LIMIT 40
-- NEW: rebuild from z-score parquet or remove submodule network visualization

-- OLD: expression
SELECT timepoint, mean_tpm, sd_tpm, n_genes FROM submodule_expression WHERE submodule=? ORDER BY ...

-- NEW: call _module_expr_by_timepoint(db, module_id)

-- OLD: GO enrichment terms
SELECT term_name FROM module_enrichment
WHERE module=? AND module_type='submodule' AND significant='TRUE' AND source IN ('GO:BP','GO:CC','GO:MF')
ORDER BY p_value LIMIT 4

-- NEW
SELECT e.term_name FROM go_module_enrichment e
JOIN module_table m ON e.module_id = m.module_id
WHERE m.module_name=? ORDER BY e.p_value LIMIT 4

-- OLD: top TFs
SELECT tf, direction FROM perturbation_edges
WHERE module=? AND module_type='submodule' ORDER BY ABS(mean_NES) DESC LIMIT 5

-- NEW
SELECT gene_name AS tf, direction FROM gsea_tf_table
WHERE module=? AND module_collection='hotspot_submodule' ORDER BY ABS(mean_NES) DESC LIMIT 5

-- OLD: neighbor supermodule lookup
SELECT id, supermodule FROM submodule_nodes WHERE id IN (...)

-- NEW: supermodule is encoded in the module_name itself (e.g. DE-3.6 → DE-3)
-- No query needed; derive with:  mod_name.rsplit('.', 1)[0]
```

---

## Gene page — `gene_page(gene_name)`

```sql
-- OLD: expression profile (gene_expression table removed)
SELECT timepoint, mean_tpm, replicates FROM gene_expression WHERE gene=? ORDER BY ...

-- NEW: two steps
--  1. resolve_gene(db, gene_name)  → dict with gene_id
--  2. _gene_expr_by_timepoint(db, gene_id)  → [{timepoint, mean_tpm, replicates}]

-- OLD: module membership
SELECT submodule, supermodule FROM submodule_genes WHERE gene=?
SELECT module FROM module_genes WHERE gene=?

-- NEW: single query returns all module memberships (sub and super)
SELECT m.module_name, m.source
FROM gene_table g
JOIN gene_module_table gm ON g.gene_id = gm.gene_id
JOIN module_table m ON gm.module_id = m.module_id
WHERE g.gene_name=?

-- OLD: TF check
SELECT 1 FROM perturbation_edges WHERE tf=?
UNION ALL SELECT 1 FROM binding_edges WHERE tf=? LIMIT 1

-- NEW
SELECT 1 FROM gsea_tf_table WHERE gene_name=? LIMIT 1
-- OR: SELECT in_perturbation_library FROM gene_table WHERE gene_name=?

-- OLD: TF description
SELECT name, summary FROM tf_descriptions WHERE tf=?

-- NEW: use gene_table fields
SELECT gene_name, Summary FROM gene_table WHERE gene_name=?
-- (tf_descriptions table removed; Summary is the gene annotation field)

-- OLD: supermodule-level perturbation (as TF)
SELECT module, mean_NES, n_sig_gRNA, padj FROM perturbation_edges
WHERE tf=? AND module_type='supermodule'

-- NEW
SELECT gene_set AS module, mean_NES, n_grnas, min_padj AS padj
FROM gsea_tf_table
WHERE gene_name=? AND module_collection='hotspot_supermodule'

-- OLD: supermodule-level binding (as TF)
SELECT module, odds_ratio FROM binding_edges
WHERE tf=? AND module_type='supermodule'

-- NEW
SELECT m.module_name AS module, MAX(e."OR") AS odds_ratio
FROM tf_module_enrichment e
JOIN module_table m ON e.module_id = m.module_id
JOIN tf_dataset_table d ON e.dataset_id = d.dataset_id
WHERE d.tf_gene_name=? AND m.source='hotspot_supermodule'
GROUP BY m.module_name

-- OLD: submodule listing for a supermodule
SELECT id, n_genes, color, within_mean_z FROM submodule_nodes WHERE supermodule=? ORDER BY id

-- NEW (note: color/within_mean_z unavailable)
SELECT module_id, module_name, size FROM module_table
WHERE source='hotspot_submodule' AND module_name LIKE ?  -- LIKE = supermod + '.%'
ORDER BY module_name

-- OLD: submodule-level perturbation (as TF)
SELECT module, mean_NES, n_sig_gRNA, padj FROM perturbation_edges
WHERE tf=? AND module_type='submodule'

-- NEW
SELECT module, mean_NES, n_grnas, min_padj AS padj
FROM gsea_tf_table
WHERE gene_name=? AND module_collection='hotspot_submodule'

-- OLD: submodule-level binding (as TF)
SELECT module, odds_ratio FROM binding_edges
WHERE tf=? AND module_type='submodule'

-- NEW
SELECT m.module_name AS module, MAX(e."OR") AS odds_ratio
FROM tf_module_enrichment e
JOIN module_table m ON e.module_id = m.module_id
JOIN tf_dataset_table d ON e.dataset_id = d.dataset_id
WHERE d.tf_gene_name=? AND m.source='hotspot_submodule'
GROUP BY m.module_name

-- OLD: submodule metadata for merged sub-ids
SELECT id, supermodule, n_genes, color FROM submodule_nodes WHERE id IN (...)

-- NEW: supermodule = mod_name.rsplit('.', 1)[0]; no color/within_mean_z
SELECT module_name, size FROM module_table
WHERE source='hotspot_submodule' AND module_name IN (...)

-- OLD: gene locus
SELECT chr, MIN(start) AS locus_start, MAX(end) AS locus_end
FROM gene_regions WHERE gene_name=? AND gene_region_subtype IN ('proximal','distal')
GROUP BY chr ORDER BY (MAX(end)-MIN(start)) DESC LIMIT 1

-- NEW (table rename only)
SELECT chr, MIN(start) AS locus_start, MAX(end) AS locus_end
FROM gene_region_table WHERE gene_name=? AND gene_region_subtype IN ('proximal','distal')
GROUP BY chr ORDER BY (MAX(end)-MIN(start)) DESC LIMIT 1
```

---

## GO page — `_go_page(term_name)`

```sql
-- OLD: resolve bare GO ID
SELECT DISTINCT term FROM go_genes WHERE term LIKE ?

-- NEW: go_term_table has the canonical name and ID separately
SELECT go_id, go_name FROM go_term_table WHERE go_id=? LIMIT 1
-- (URL scheme should move to go_id-based; old term-name-as-key scheme is fragile)

-- OLD: member genes
SELECT gene FROM go_genes WHERE term=? ORDER BY gene

-- NEW
SELECT g.gene_name FROM go_gene_table gg
JOIN gene_table g USING(gene_id)
JOIN go_term_table gt ON gg.go_id = gt.go_id
WHERE gt.go_name=? OR gt.go_id=?
ORDER BY g.gene_name

-- OLD: enriched modules
SELECT module, module_type, p_value, term_size, intersection_size,
       precision_val, recall, query_size
FROM module_enrichment WHERE term_id=? AND source='GO:BP' ORDER BY p_value

-- NEW
SELECT m.module_name, m.source, e.p_value, e.term_size, e.intersection_size,
       e.precision, e.recall, e.query_size
FROM go_module_enrichment e
JOIN module_table m ON e.module_id = m.module_id
WHERE e.go_id=? ORDER BY e.p_value

-- OLD: TFs perturbing this GO term
SELECT tf, mean_NES, n_sig_gRNA, padj FROM perturbation_edges
WHERE module=? AND module_type='go' ORDER BY ABS(mean_NES) DESC

-- NEW: gsea_tf_table has go_id FK; check gene_set_collection for GO entries
SELECT gene_name AS tf, mean_NES, n_grnas, min_padj AS padj
FROM gsea_tf_table
WHERE go_id=? ORDER BY ABS(mean_NES) DESC
-- NOTE: verify that GO terms are populated in gsea_tf_table

-- OLD: TFs binding this GO term gene set
SELECT tf, odds_ratio FROM binding_edges WHERE module=? AND module_type='go'

-- NEW: gmt_enrichment covers TF target enrichment vs gene sets (may include GO sets)
SELECT d.tf_gene_name AS tf, MAX(e."OR") AS odds_ratio
FROM gmt_enrichment e
JOIN tf_dataset_table d ON e.dataset_id = d.dataset_id
WHERE e.gene_set=? AND e.gene_set_collection LIKE '%GO%'
GROUP BY d.tf_gene_name
-- NOTE: verify gene_set_collection values for GO entries in gmt_enrichment
```

---

## Element page — `element_page(atac_peak_id)`

```sql
-- OLD: peak lookup
SELECT atac_peak_id, chr, start, end, atac_peak_name FROM atac_peaks WHERE atac_peak_id=?

-- NEW (table rename)
SELECT atac_peak_id, chr, start, end, atac_peak_name FROM atac_peak_table WHERE atac_peak_id=?

-- OLD: nearby peaks
SELECT chr, start, end, atac_peak_name FROM atac_peaks
WHERE chr=? AND start>=? AND end<=? AND atac_peak_id!=? ORDER BY start

-- NEW
SELECT chr, start, end, atac_peak_name FROM atac_peak_table
WHERE chr=? AND start>=? AND end<=? AND atac_peak_id!=? ORDER BY start

-- Query helpers _query_element_tfs, _query_element_genes, _query_tf_gene_link:
-- All tf_datasets  → tf_dataset_table
-- All atac_peaks   → atac_peak_table
-- All gene_regions → gene_region_table
-- (No structural changes beyond table renames in these joins)

-- NEW: also query multiome_atac_overlaps for gene-peak links
SELECT g.gene_name, o.cell_type, o.link_type, o.distance_to_tss
FROM multiome_atac_overlaps o
JOIN gene_table g ON o.gene_id = g.gene_id
WHERE o.atac_peak_id=?
```

---

## TF-gene link page — `tf_gene_link_page(tf, gene)`

```sql
-- OLD: expression (both tf and gene)
SELECT timepoint, mean_tpm FROM gene_expression WHERE gene=? ORDER BY ...

-- NEW: resolve_gene(db, name) → _gene_expr_by_timepoint(db, gene_id)

-- OLD: TF description display name
SELECT name FROM tf_descriptions WHERE tf=?

-- NEW: gene_table.gene_name is already the display name; no separate table
-- Use: SELECT gene_name FROM gene_table WHERE gene_name=?
-- OR: gene_table.Summary for a description string

-- NEW: tf_gene_links provides direct evidence (use alongside _query_tf_gene_link):
SELECT tgl.source, tgl.gene_name
FROM tf_gene_links tgl
JOIN tf_dataset_table d ON tgl.dataset_id = d.dataset_id
WHERE d.tf_gene_name=? AND tgl.gene_name=?
```

---

## Dataset page — `dataset_page(dataset_id)`

```sql
-- OLD: dataset lookup
SELECT * FROM tf_datasets WHERE dataset_id=?

-- NEW
SELECT * FROM tf_dataset_table WHERE dataset_id=?

-- OLD: TF description
SELECT name, summary FROM tf_descriptions WHERE tf=?

-- NEW: gene_table.Summary
SELECT gene_name, Summary FROM gene_table WHERE gene_name=?

-- OLD: peak count (unchanged)
SELECT COUNT(*) FROM tf_peaks WHERE dataset_id=?

-- OLD: elements overlapped (unchanged join)
SELECT COUNT(DISTINCT ato.atac_peak_id)
FROM atac_tf_overlaps ato JOIN tf_peaks tp ON ato.peak_id=tp.peak_id
WHERE tp.dataset_id=?

-- OLD: target gene count (expensive join)
SELECT COUNT(DISTINCT gr.gene_name)
FROM tf_gene_overlaps tgo JOIN tf_peaks tp ON tgo.peak_id=tp.peak_id
JOIN gene_regions gr ON tgo.gene_region_id=gr.gene_region_id
WHERE tp.dataset_id=?

-- NEW: use tf_dataset_gene_table (pre-computed, much cheaper)
SELECT COUNT(DISTINCT gene_name) FROM tf_dataset_gene_table WHERE dataset_id=?

-- OLD: gene rows (expensive join)
SELECT gr.gene_name, gr.gene_region_subtype, COUNT(DISTINCT tp.peak_id) AS peak_count
FROM tf_gene_overlaps tgo JOIN tf_peaks tp ON tgo.peak_id=tp.peak_id
JOIN gene_regions gr ON tgo.gene_region_id=gr.gene_region_id
WHERE tp.dataset_id=?
GROUP BY gr.gene_name, gr.gene_region_subtype

-- NEW: use tf_dataset_gene_table (simpler)
SELECT gene_name, match_type FROM tf_dataset_gene_table WHERE dataset_id=?

-- OLD: related datasets
SELECT dataset_id, dataset, cell_type, cell_type_group, source
FROM tf_datasets WHERE tf_gene_name=? AND dataset_id!=? ORDER BY cell_type, dataset

-- NEW
SELECT dataset_id, dataset, cell_type, cell_type_group, source
FROM tf_dataset_table WHERE tf_gene_name=? AND dataset_id!=? ORDER BY cell_type, dataset

-- OLD: element rows
SELECT ap.atac_peak_id, ap.chr, ap.start, ap.end, ...
FROM atac_tf_overlaps ato JOIN tf_peaks tp ... JOIN atac_peaks ap ... JOIN tf_gene_overlaps ...
WHERE tp.dataset_id=?

-- NEW: atac_peaks → atac_peak_table  (rest unchanged)
```

---

## API: `/api/edges`

```sql
-- OLD: perturbation edges (all)
SELECT tf, module, mean_NES, n_sig_gRNA, padj FROM perturbation_edges

-- NEW
SELECT gene_name AS tf, gene_set AS module, mean_NES, n_grnas, min_padj AS padj,
       direction, module_collection
FROM gsea_tf_table
WHERE module_collection IN ('hotspot_supermodule', 'hotspot_submodule')

-- OLD: binding edges (all)
SELECT tf, module, odds_ratio, padj FROM binding_edges

-- NEW
SELECT d.tf_gene_name AS tf, m.module_name AS module,
       MAX(e."OR") AS odds_ratio, MIN(e.padj_fisher) AS padj
FROM tf_module_enrichment e
JOIN tf_dataset_table d ON e.dataset_id = d.dataset_id
JOIN module_table m ON e.module_id = m.module_id
GROUP BY d.tf_gene_name, m.module_name
```

---

## API: `/api/tf/<tf_name>`

```sql
-- OLD: perturbation
SELECT module, mean_NES, n_sig_gRNA, padj FROM perturbation_edges WHERE tf=?

-- NEW
SELECT gene_set AS module, module_collection, mean_NES, n_grnas, min_padj AS padj, direction
FROM gsea_tf_table WHERE gene_name=?

-- OLD: binding
SELECT module, odds_ratio, log2or, padj FROM binding_edges WHERE tf=?

-- NEW
SELECT m.module_name AS module, MAX(e."OR") AS odds_ratio, MIN(e.padj_fisher) AS padj,
       m.source AS module_collection
FROM tf_module_enrichment e
JOIN tf_dataset_table d ON e.dataset_id = d.dataset_id
JOIN module_table m ON e.module_id = m.module_id
WHERE d.tf_gene_name=?
GROUP BY m.module_name
```

---

## API: `/api/module/<mod_name>`

```sql
-- OLD: perturbation
SELECT tf, mean_NES, n_sig_gRNA, padj FROM perturbation_edges WHERE module=?

-- NEW
SELECT gene_name AS tf, mean_NES, n_grnas, min_padj AS padj, direction
FROM gsea_tf_table WHERE module=?

-- OLD: binding
SELECT tf, odds_ratio, log2or, padj FROM binding_edges WHERE module=?

-- NEW
SELECT d.tf_gene_name AS tf, MAX(e."OR") AS odds_ratio, MIN(e.padj_fisher) AS padj
FROM tf_module_enrichment e
JOIN tf_dataset_table d ON e.dataset_id = d.dataset_id
JOIN module_table m ON e.module_id = m.module_id
WHERE m.module_name=?
GROUP BY d.tf_gene_name
```

---

## API: `/api/module/<mod_name>/detail`

```sql
-- OLD: GO enrichment
SELECT term_id, source, term_name, p_value, term_size, intersection_size, precision_val, recall
FROM module_enrichment WHERE module=? AND source IN ('GO:BP','GO:CC','GO:MF')

-- NEW
SELECT e.go_id AS term_id, e.source, e.term_name, e.p_value,
       e.term_size, e.intersection_size, e.precision, e.recall
FROM go_module_enrichment e
JOIN module_table m ON e.module_id = m.module_id
WHERE m.module_name=? ORDER BY e.p_value

-- OLD: gene list
SELECT gene FROM module_genes WHERE module=? ORDER BY gene

-- NEW
SELECT gene_name FROM gene_module_table
JOIN module_table USING(module_id)
WHERE module_name=? AND gene_name IS NOT NULL ORDER BY gene_name

-- OLD: submodule colors
SELECT id, color FROM submodule_nodes WHERE supermodule=?

-- NEW: color not in new schema; derive from MODULE_COLORS using parent name

-- OLD: gene→submodule mapping
SELECT gene, submodule FROM submodule_genes WHERE submodule IN (SELECT id FROM submodule_nodes WHERE supermodule=?)

-- NEW
SELECT gm_sub.gene_name, m_sub.module_name AS submodule
FROM module_table m_super
JOIN module_table m_sub ON m_sub.source='hotspot_submodule' AND m_sub.module_name LIKE m_super.module_name || '.%'
JOIN gene_module_table gm_sub ON gm_sub.module_id = m_sub.module_id
WHERE m_super.module_name=? AND gm_sub.gene_name IS NOT NULL
```

---

## API: `/api/module/<mod_name>/enrichment`

```sql
-- OLD
SELECT term_id, term_name, source, p_value, term_size, intersection_size, precision_val, recall
FROM module_enrichment WHERE module=? AND source IN ('GO:BP','GO:CC','GO:MF') ORDER BY p_value

-- NEW
SELECT e.go_id AS term_id, e.term_name, e.source, e.p_value,
       e.term_size, e.intersection_size, e.precision, e.recall
FROM go_module_enrichment e
JOIN module_table m ON e.module_id = m.module_id
WHERE m.module_name=? ORDER BY e.p_value
```

---

## API: `/api/search`

```sql
-- OLD: active TFs
SELECT DISTINCT tf FROM tf_descriptions WHERE tf LIKE ? ORDER BY tf LIMIT 8

-- NEW
SELECT gene_name FROM gene_table
WHERE gene_name LIKE ? AND in_perturbation_library=1 ORDER BY gene_name LIMIT 8

-- OLD: modules
SELECT DISTINCT module FROM module_descriptions WHERE module LIKE ? ORDER BY module LIMIT 5

-- NEW
SELECT module_name FROM module_table
WHERE module_name LIKE ? AND source='hotspot_supermodule' AND module_name != 'unassigned'
ORDER BY module_name LIMIT 5

-- OLD: submodules
SELECT DISTINCT id FROM submodule_nodes WHERE id LIKE ? ORDER BY id LIMIT 5

-- NEW
SELECT module_name FROM module_table
WHERE module_name LIKE ? AND source='hotspot_submodule'
ORDER BY module_name LIMIT 5

-- OLD: GO terms
SELECT DISTINCT term FROM go_genes WHERE term LIKE ? ORDER BY LENGTH(term) LIMIT 5

-- NEW (search by name; return go_id for URL routing)
SELECT go_id, go_name FROM go_term_table
WHERE go_name LIKE ? AND is_obsolete=0 ORDER BY LENGTH(go_name) LIMIT 5
-- NOTE: GO page URLs should switch from term-name to go_id

-- OLD: non-TF genes
SELECT DISTINCT gene FROM gene_expression WHERE gene LIKE ? ORDER BY gene LIMIT 20

-- NEW
SELECT gene_name FROM gene_table WHERE gene_name LIKE ? ORDER BY gene_name LIMIT 20

-- NEW: synonym search (not available before)
SELECT DISTINCT g.gene_name FROM gene_synonym s
JOIN gene_table g USING(gene_id)
WHERE s.synonym LIKE ? ORDER BY g.gene_name LIMIT 10
```

---

## API: `/api/link-tfs`

```sql
-- OLD
SELECT DISTINCT tf_gene_name FROM tf_datasets WHERE tf_gene_name IS NOT NULL ORDER BY tf_gene_name

-- NEW (table rename only)
SELECT DISTINCT tf_gene_name FROM tf_dataset_table WHERE tf_gene_name IS NOT NULL ORDER BY tf_gene_name
```

---

## API: `/api/link-genes`

```sql
-- OLD
SELECT DISTINCT gene_name FROM gene_regions WHERE gene_name IS NOT NULL ORDER BY gene_name

-- NEW (table rename only)
SELECT DISTINCT gene_name FROM gene_region_table WHERE gene_name IS NOT NULL ORDER BY gene_name
```

---

## API: `/api/submodule/<name>` and `/api/submodule/<name>/gene_network`

```sql
-- OLD: submodule node
SELECT * FROM submodule_nodes WHERE id=?

-- NEW
SELECT module_id, module_name, source, size FROM module_table
WHERE module_name=? AND source='hotspot_submodule'

-- OLD: genes
SELECT gene FROM submodule_genes WHERE submodule=? ORDER BY gene

-- NEW
SELECT gene_name FROM gene_module_table WHERE module_id=? AND gene_name IS NOT NULL ORDER BY gene_name

-- OLD: neighbors (submodule edges — REMOVED from new schema)
SELECT source, target, mean_z, n_pairs FROM submodule_edges
WHERE source=? OR target=? ORDER BY ABS(mean_z) DESC LIMIT 40

-- NEW: no equivalent; must rebuild from z-score parquet or remove endpoint

-- OLD: gene-gene edges (submodule_gene_edges — NOT in TfDatabase schema)
SELECT gene1, gene2, z_score FROM submodule_gene_edges WHERE submodule=? ORDER BY z_score DESC

-- NEW: rebuild from parquet at startup and store in-process, or add to new DB
```

---

## API: `/api/module/<mod_name>/submodule_network`

```sql
-- OLD
SELECT id, n_genes, color, within_mean_z FROM submodule_nodes WHERE supermodule=?
SELECT source, target, mean_z, n_pairs FROM submodule_edges WHERE source LIKE ? OR target LIKE ?

-- NEW: submodule_edges table removed; network data must be rebuilt
-- Nodes (partial — no color/within_mean_z):
SELECT module_id, module_name, size FROM module_table
WHERE source='hotspot_submodule' AND module_name LIKE ?  -- LIKE = mod_name + '.%'
```

---

## Tables renamed only (no structural changes)

| Old name | New name |
|---|---|
| `tf_datasets` | `tf_dataset_table` |
| `atac_peaks` | `atac_peak_table` |
| `gene_regions` | `gene_region_table` |

All joins using these tables need their FROM/JOIN clauses updated. Column names are unchanged.

---

## Tables removed (no equivalent in TfDatabase)

| Old table | Notes |
|---|---|
| `perturbation_edges` | → `gsea_tf_table` (more detailed) |
| `binding_edges` | → `tf_module_enrichment` (more detailed) |
| `module_genes` | → `gene_module_table` (generic) |
| `submodule_nodes` | → `module_table WHERE source='hotspot_submodule'` (no color/within_mean_z) |
| `submodule_genes` | → `gene_module_table` |
| `submodule_edges` | → **no equivalent** — rebuild or remove |
| `submodule_gene_edges` | → **no equivalent** — rebuild or remove |
| `module_enrichment` | → `go_module_enrichment` (uses module_id FK, not name) |
| `go_genes` | → `go_gene_table` (uses gene_id FK) |
| `module_descriptions` | → **no equivalent** — keep TSV externally |
| `tf_descriptions` | → `gene_table.Summary` |
| `gene_expression` | → `bulk_expression` (per-sample; use `_gene_expr_by_timepoint()`) |
| `module_expression` | → computed via `_module_expr_by_timepoint()` |
| `submodule_expression` | → computed via `_module_expr_by_timepoint()` |
| `expressed_tfs` | → `gene_table.in_perturbation_library` or `gsea_tf_table` |

## Tables new in TfDatabase (not in old schema)

`gene_table`, `gene_tss_table`, `gene_synonym`, `tf_dataset_gene_table`, `tf_gene_links`,
`multiome_linked_peaks`, `multiome_atac_overlaps`, `grna_table`, `de_results`,
`module_table`, `gene_module_table`, `go_term_table`, `go_gene_table`,
`tf_module_enrichment`, `gmt_enrichment`, `gsea_grna_table`, `gsea_tf_table`,
`bulk_expression`, `go_module_enrichment`
