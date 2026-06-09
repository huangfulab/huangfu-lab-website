# Flask Refactor Plan

**Date:** 2026-06-09 08:30:06
**Mode:** Auto
**Files analyzed:** 2 (`app.py`, `perturbseq_bp.py`)
**Redundancy groups found:** 7 (3 bug-risk, 2 maintenance, 2 style/perf)

---

## Summary

The app is a single-blueprint Flask app (`perturbseq_bp.py`, 3414 lines) backed by a single SQLite database. It is mostly clean, but contains three critical redundancy bugs: a duplicate timepoint-order constant that will silently diverge if ever updated, a ghost column name in a sort map that crashes the server when used, and a dead `if/else` with identical branches that hides missing logic. Two maintenance-level duplications and two style issues round out the findings.

---

## Redundancy Groups

### Group 1: `_SPARK_TP_ORDER` duplicates `_TP_ORDER` exactly

**Type:** Copy-paste drift
**Risk:** 🔴 Bug risk (copies will eventually drift)
**Files affected:** `perturbseq_bp.py`

**Description:**
Two constants are defined with the same value at line 65 and line 380:

```python
# line 65
_TP_ORDER = ['ES_0h', 'DE_12h', 'DE_24h', 'DE_36h', 'DE_48h', 'DE_60h', 'DE_72h']

# line 380 — IDENTICAL
_SPARK_TP_ORDER = ['ES_0h', 'DE_12h', 'DE_24h', 'DE_36h', 'DE_48h', 'DE_60h', 'DE_72h']
```

`_TP_ORDER` drives DB queries and profile computations; `_SPARK_TP_ORDER` drives sparkline rendering. If a new timepoint is added to `_TP_ORDER` (e.g., `'DE_96h'`), `_SPARK_TP_ORDER` will silently be left behind, causing sparklines to render without the new timepoint while everything else shows it.

**Proposed fix:**
Replace the definition of `_SPARK_TP_ORDER` with an alias: `_SPARK_TP_ORDER = _TP_ORDER`. The constant stays in place so call sites require no changes.

**Files to modify:** `perturbseq_bp.py` line 380

---

### Group 2: Ghost `n_ctg` column in sort map causes SQL error

**Type:** Copy-paste drift (stale rename)
**Risk:** 🔴 Bug risk (SQL error on sort)
**Files affected:** `perturbseq_bp.py`

**Description:**
`_query_tf_linked_genes_paged` has an `order_map` at line 1607–1615:

```python
order_map = {
    0: f'gt.gene_name {safe_dir}',
    1: f'sub.module_name {safe_dir}',
    2: f'n_datasets {safe_dir}',
    3: f'n_ctg {safe_dir}',          # ← column does not exist in the SELECT
    4: f'gdp.signed_pct_rank {safe_dir}',
    5: f'gdp.mean_coef {safe_dir}',
}
```

The SELECT statement aliases the same concept as `n_linked_peaks` (line 1742: `COALESCE(pc.n_linked_peaks, 0) AS n_linked_peaks`). The sort map was never updated after the rename. If the frontend sorts by column index 3, SQLite raises `no such column: n_ctg`.

**Proposed fix:**
Change `3: f'n_ctg {safe_dir}'` → `3: f'n_linked_peaks {safe_dir}'`.

**Files to modify:** `perturbseq_bp.py` line 1611

---

### Group 3: Dead `if/else` with identical branches in `gene_page` tooltip builder

**Type:** Copy-paste drift
**Risk:** 🔴 Bug risk (dead code masking missing feature)
**Files affected:** `perturbseq_bp.py`

**Description:**
Inside `gene_page` (lines 1395–1404), a `for` loop builds `module_tooltips` with a branch that checks `_r['source']` but both branches produce the same dict:

```python
# lines 1395–1404
if _r['source'] == 'hotspot_submodule':
    module_tooltips[_r['module_name']] = {
        'title': _r['title'],
        'desc': _r['standard'],
    }
else:
    module_tooltips[_r['module_name']] = {
        'title': _r['title'],
        'desc': _r['standard'],
    }
```

Compare to `_super_page` (lines 865–886) where a similar loop adds `spark_line`/`spark_area` to the tooltip for submodules in the `else` branch. The `gene_page` copy was written intending the same treatment but the extra logic was never added — leaving identical branches that mislead any future developer into thinking the distinction is intentional.

**Proposed fix:**
Collapse the `if/else` into a single unconditional assignment:
```python
module_tooltips[_r['module_name']] = {
    'title': _r['title'],
    'desc': _r['standard'],
}
```
Also drop `mt.source` from the SELECT since it is no longer used.

**Files to modify:** `perturbseq_bp.py` lines 1385–1405

---

### Group 4: Duplicated `sources`/`datasets_str` parsing in two API functions

**Type:** Repeated inline pattern
**Risk:** 🟡 Maintenance burden
**Files affected:** `perturbseq_bp.py`

**Description:**
The same 8-line block appears identically in two places:

```python
# _query_tf_linked_genes_paged, lines 1760–1768
sources = sorted(set((r['sources_str'] or '').split(','))) if r['sources_str'] else []
seen_ds = set()
datasets = []
for entry in (r['datasets_str'] or '').split(','):
    parts = entry.split('|||')
    if len(parts) == 4 and parts[0] not in seen_ds:
        seen_ds.add(parts[0])
        datasets.append({'id': parts[0], 'name': parts[1], 'cell_type': parts[2], 'source': parts[3]})

# api_tf_linked_genes_all, lines 1895–1903 — IDENTICAL
sources = sorted(set((r['sources_str'] or '').split(','))) if r['sources_str'] else []
seen_ds = set()
datasets = []
for entry in (r['datasets_str'] or '').split(','):
    parts = entry.split('|||')
    if len(parts) == 4 and parts[0] not in seen_ds:
        seen_ds.add(parts[0])
        datasets.append({'id': parts[0], 'name': parts[1], 'cell_type': parts[2], 'source': parts[3]})
```

Any change to the `|||` encoding format or the output dict shape must be applied in two places.

**Proposed fix:**
Extract a private helper `_parse_sources_datasets(r) -> tuple[list, list]` placed near `rows_to_dicts` (after line 377). Both call sites are replaced with `sources, datasets = _parse_sources_datasets(r)`.

**Files to modify:** `perturbseq_bp.py` lines ~377, 1760–1768, 1895–1903

---

### Group 5: `BW_EXTRA_DIR` and `BW_DIR` are set to the same path

**Type:** Copy-paste drift
**Risk:** 🟡 Maintenance burden
**Files affected:** `perturbseq_bp.py`

**Description:**
Lines 45–46 define two directory constants pointing to the same path:

```python
BW_DIR       = str(Path(__file__).resolve().parent / "data" / "bw" / "atac")
BW_EXTRA_DIR = str(Path(__file__).resolve().parent / "data" / "bw" / "atac")
```

`serve_bw` and `serve_bw_extra` (lines 3055 and 3061) serve files from the same folder. If the extra-dir path is ever changed, the copy must be updated too.

**Proposed fix:**
Replace the `BW_EXTRA_DIR` assignment with `BW_EXTRA_DIR = BW_DIR`.

**Files to modify:** `perturbseq_bp.py` line 46

---

### Group 6: Heavy imports inside `_coef_kde_data_db` run on every call

**Type:** Import hygiene
**Risk:** 🔵 Style/perf
**Files affected:** `perturbseq_bp.py`

**Description:**
Lines 2968–2969 import `numpy` and `scipy.stats.gaussian_kde` inside a route-backing function:

```python
def _coef_kde_data_db(db, tf: str, gene: str, n_pts: int = 100):
    import numpy as np
    from scipy.stats import gaussian_kde
```

Python caches imports after the first call, so the performance hit is minor, but the convention breaks IDE analysis, type checkers, and linters.

**Proposed fix:**
Move both imports to the module top-level block (after line 9).

**Files to modify:** `perturbseq_bp.py` lines 2968–2969 and top-of-file

---

### Group 7: `get_tc_db` is a dead alias for `get_db`

**Type:** Duplicate function
**Risk:** 🔵 Style/perf
**Files affected:** `perturbseq_bp.py`

**Description:**
Lines 303–305 define a wrapper that only forwards to `get_db`:

```python
def get_tc_db():
    return get_db()
```

`get_tc_db` is not imported by `app.py` or any other file and is never called within `perturbseq_bp.py` itself. It is dead code.

**Proposed fix:**
Remove the `get_tc_db` function entirely.

**Files to modify:** `perturbseq_bp.py` lines 303–305
