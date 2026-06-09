# Flask Refactor Fix Report

**Date:** 2026-06-09 10:07:50
**Mode:** Auto
**Plan source:** `logs/flask_refactor_plan_2026-06-09_100433.md`
**Files modified:** `perturbseq_bp.py`
**Groups fixed:** 5 of 5 (0 bug-risk, 3 maintenance, 2 style/perf)

---

## Summary

All five structural redundancies from the second-pass plan were resolved. The most impactful changes extract a `_query_db_raw` helper (eliminating three separate raw SQLite connection boilerplate blocks), consolidate four inline gene-sorting patterns into the existing `_top_genes_by_pub` helper, and extract a `_attach_sparklines` helper that replaces a tripled loop. Two dead functions (`_load_gene_name_map` and `_tp_order`) were removed entirely.

---

## Fixes Applied

### Group 1: Tripled `sqlite3.connect / execute / close` boilerplate in cache loaders — 🟡 Maintenance burden

**What was wrong:**
`_load_gene_pub_counts`, `_load_expressed_tfs`, and `_load_perturbed_tfs` each opened a raw SQLite connection with an identical try/finally structure. Any change to connection setup (e.g., adding a pragma, timeout, or row_factory) required three synchronised edits.

**What changed:**

```python
# Before — each loader had this structure independently:
conn = sqlite3.connect(DB_PATH)
try:
    rows = conn.execute("SELECT ...").fetchall()
finally:
    conn.close()

# After — new helper placed before _load_lambert_tfs:
def _query_db_raw(sql: str, params=()) -> list:
    conn = sqlite3.connect(DB_PATH)
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()

# Each loader now uses:
rows = _query_db_raw("SELECT ...")
```

**Why this matters:**
SQLite connection setup is now maintained in one place; all three cache loaders stay in sync automatically.

---

### Group 2: `_top_genes_by_pub` helper bypassed by 4 identical inline sorts — 🟡 Maintenance burden

**What was wrong:**
Four call sites in `all_modules_page` (lines ~471, ~517, ~578) and `_super_page` (line ~815) sorted genes inline with an identical lambda instead of using the existing `_top_genes_by_pub` helper. Any change to the sort key (e.g., secondary tie-break) required five synchronised edits.

**What changed:**

```python
# Before (×4):
notable = sorted(genes, key=lambda g: pub_counts.get(g, 0), reverse=True)[:5]
# and:
s['top_genes'] = sorted(genes, key=lambda g: pub_counts.get(g, 0), reverse=True)[:5]

# After (×4):
notable = _top_genes_by_pub(genes, n=5)
# and:
s['top_genes'] = _top_genes_by_pub(genes, n=5)
```

**Why this matters:**
The sort logic now lives exclusively in `_top_genes_by_pub`; all four call sites track it automatically.

---

### Group 3: Tripled sparkline-attachment loop in `all_modules_page` — 🟡 Maintenance burden

**What was wrong:**
Three structurally identical loops in `all_modules_page` (lines ~600–617) each iterated a list of module dicts, looked up the expression map, and called `_spark_pts`. The only variation was the dict key used for the map lookup (`'id'` vs `'name'`).

**What changed:**

```python
# Before (×3 with minor variation):
for c in clusters:
    tp_dict = gc_expr_map.get(c['id'], {})
    vals = [tp_dict.get(tp) for tp in _SPARK_TP_ORDER]
    c['spark_line'], c['spark_area'] = _spark_pts(vals)

# New helper added after _spark_pts:
def _attach_sparklines(items, expr_map, key='name'):
    for item in items:
        tp_dict = expr_map.get(item[key], {})
        vals = [tp_dict.get(tp) for tp in _SPARK_TP_ORDER]
        item['spark_line'], item['spark_area'] = _spark_pts(vals)

# After — three calls replace three loops:
_attach_sparklines(clusters, gc_expr_map, key='id')
_attach_sparklines(supermodules, sup_expr_map)
_attach_sparklines(submodules, sub_expr_map)
```

**Why this matters:**
Any change to sparkline attachment (e.g., adding a third SVG path or switching to `_SPARK_TP_ORDER` indices) is now made in one place.

---

### Group 4: Dead function `_load_gene_name_map` — 🔵 Style/perf

**What was wrong:**
`_load_gene_name_map` loaded a TSV into a `_gene_name_map` module-level cache but was never called anywhere in the file. The accompanying `_gene_name_map: dict | None = None` global was also dead.

**What changed:**
- Removed `_gene_name_map: dict | None = None` from module-level globals
- Removed `_load_gene_name_map` function entirely

**Why this matters:**
Dead code is gone; developers no longer need to wonder whether the TSV-backed name map is relevant to the live codebase.

---

### Group 5: Dead function `_tp_order` — 🔵 Style/perf

**What was wrong:**
`_tp_order()` returned a hardcoded SQL `ORDER BY CASE timepoint …` string but was never called anywhere in the file. The same ordering is handled elsewhere by `_TP_ORDER` and direct SQL.

**What changed:**
Removed `_tp_order` function entirely.

**Why this matters:**
Removes a misleading function that suggested a centralised timepoint-ordering path that did not actually exist in any active query.

---

## Groups Not Fixed

None — all groups from the plan were resolved.

---

## File Change Index

| File | Nature of change |
|------|-----------------|
| `perturbseq_bp.py` | Added `_query_db_raw(sql, params)` helper before `_load_lambert_tfs` |
| `perturbseq_bp.py` | `_load_gene_pub_counts`, `_load_expressed_tfs`, `_load_perturbed_tfs` refactored to use `_query_db_raw` |
| `perturbseq_bp.py` | Removed `_load_gene_name_map` function and `_gene_name_map` global |
| `perturbseq_bp.py` | 3 inline sorts in `all_modules_page` replaced with `_top_genes_by_pub(genes, n=5)` |
| `perturbseq_bp.py` | 1 inline sort in `_super_page` replaced with `_top_genes_by_pub(genes, n=5)` |
| `perturbseq_bp.py` | Added `_attach_sparklines(items, expr_map, key='name')` helper after `_spark_pts` |
| `perturbseq_bp.py` | 3 sparkline loops in `all_modules_page` replaced with `_attach_sparklines(...)` calls |
| `perturbseq_bp.py` | Removed dead `_tp_order()` function |
