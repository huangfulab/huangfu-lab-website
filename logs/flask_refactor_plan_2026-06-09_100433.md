# Flask Refactor Plan

**Date:** 2026-06-09 10:04:33
**Mode:** Auto
**Files analyzed:** 2 (`app.py`, `perturbseq_bp.py`)
**Redundancy groups found:** 5 (0 bug-risk, 0 likely-bug, 3 maintenance, 2 style/perf)

---

## Summary

The previous pass eliminated all correctness bugs. This pass targets structural redundancies that remain: a tripled raw-DB-connection boilerplate in the three module-level cache loaders, four call sites that bypass the existing `_top_genes_by_pub` helper with identical inline sorts, a tripled sparkline-attachment loop in `all_modules_page`, and two dead functions that were never called. No new correctness issues were found.

---

## Redundancy Groups

### Group 1: Tripled `sqlite3.connect / execute / close` boilerplate in cache loaders

**Type:** Repeated inline pattern
**Risk:** 🟡 Maintenance burden
**Files affected:** `perturbseq_bp.py`

**Description:**
`_load_gene_pub_counts` (lines 271–287), `_load_expressed_tfs` (lines 307–322), and `_load_perturbed_tfs` (lines 325–341) all open a raw SQLite connection with an identical try/finally structure:

```python
# _load_expressed_tfs (lines 310–320)
conn = sqlite3.connect(DB_PATH)
try:
    rows = conn.execute(
        "SELECT gene_name FROM gene_table WHERE in_perturbation_library=1"
    ).fetchall()
finally:
    conn.close()

# _load_perturbed_tfs (lines 328–338) — IDENTICAL STRUCTURE
conn = sqlite3.connect(DB_PATH)
try:
    rows = conn.execute(
        "SELECT DISTINCT gene_name FROM gsea_tf_table ..."
    ).fetchall()
finally:
    conn.close()
```

If the connection setup ever needs to change (e.g., adding `row_factory`, pragma, or timeout), it must be updated in three places.

**Proposed fix:**
Extract a private `_query_db_raw(sql, params=())` helper that opens a connection, executes once, closes, and returns the rows. Place it immediately before `_load_lambert_tfs`. Replace all three try/connect/close blocks with a call to it.

**Files to modify:** `perturbseq_bp.py` lines ~257–341

---

### Group 2: `_top_genes_by_pub` helper bypassed by 4 identical inline sorts

**Type:** Repeated inline pattern (helper exists but not used)
**Risk:** 🟡 Maintenance burden
**Files affected:** `perturbseq_bp.py`

**Description:**
`_top_genes_by_pub(genes, n=8)` exists at line 70 but four call sites in `all_modules_page` and `_super_page` sort inline instead:

```python
# all_modules_page, line 493 (clusters)
notable = sorted(genes, key=lambda g: pub_counts.get(g, 0), reverse=True)[:5]

# all_modules_page, line 539 (supermodules) — IDENTICAL
notable = sorted(genes, key=lambda g: pub_counts.get(g, 0), reverse=True)[:5]

# all_modules_page, line 600 (submodules) — IDENTICAL
notable = sorted(genes, key=lambda g: pub_counts.get(g, 0), reverse=True)[:5]

# _super_page, line 837 (submodule top_genes) — IDENTICAL
s['top_genes'] = sorted(genes, key=lambda g: pub_counts.get(g, 0), reverse=True)[:5]
```

Any change to the sort key (e.g., secondary tie-break) must be applied in five places.

**Proposed fix:**
Replace all four inline sorts with `_top_genes_by_pub(genes, n=5)`. The helper calls `_load_gene_pub_counts()` internally, which returns the same cached global dict as the local `pub_counts` variable.

**Files to modify:** `perturbseq_bp.py` lines 493, 539, 600, 837

---

### Group 3: Tripled sparkline-attachment loop in `all_modules_page`

**Type:** Repeated inline pattern
**Risk:** 🟡 Maintenance burden
**Files affected:** `perturbseq_bp.py`

**Description:**
`all_modules_page` (lines 622–639) contains three structurally identical loops that attach `spark_line` / `spark_area` to module dicts:

```python
# clusters (lines 622–625)
for c in clusters:
    tp_dict = gc_expr_map.get(c['id'], {})
    vals = [tp_dict.get(tp) for tp in _SPARK_TP_ORDER]
    c['spark_line'], c['spark_area'] = _spark_pts(vals)

# supermodules (lines 628–632) — IDENTICAL except key='name'
for m in supermodules:
    tp_dict = sup_expr_map.get(m['name'], {})
    vals = [tp_dict.get(tp) for tp in _SPARK_TP_ORDER]
    m['spark_line'], m['spark_area'] = _spark_pts(vals)

# submodules (lines 635–639) — IDENTICAL
for m in submodules:
    tp_dict = sub_expr_map.get(m['name'], {})
    vals = [tp_dict.get(tp) for tp in _SPARK_TP_ORDER]
    m['spark_line'], m['spark_area'] = _spark_pts(vals)
```

**Proposed fix:**
Extract `_attach_sparklines(items, expr_map, key='name')` placed near `_spark_pts`. Replace the three loops with:
```python
_attach_sparklines(clusters, gc_expr_map, key='id')
_attach_sparklines(supermodules, sup_expr_map)
_attach_sparklines(submodules, sub_expr_map)
```

**Files to modify:** `perturbseq_bp.py` lines ~424, 622–639

---

### Group 4: Dead function `_load_gene_name_map`

**Type:** Duplicate function (dead code)
**Risk:** 🔵 Style/perf
**Files affected:** `perturbseq_bp.py`

**Description:**
`_load_gene_name_map` (lines 290–303) loads a TSV file into a module-level cache but is never called anywhere in the file or imported by any other file.

**Proposed fix:**
Remove the function and the `_gene_name_map: dict | None = None` global at line 254.

**Files to modify:** `perturbseq_bp.py` lines 254, 290–303

---

### Group 5: Dead function `_tp_order()`

**Type:** Duplicate function (dead code)
**Risk:** 🔵 Style/perf
**Files affected:** `perturbseq_bp.py`

**Description:**
`_tp_order()` (lines 689–693) returns a hardcoded SQL `ORDER BY CASE timepoint …` clause string, but is never called anywhere.

**Proposed fix:**
Remove the function entirely.

**Files to modify:** `perturbseq_bp.py` lines 689–693
