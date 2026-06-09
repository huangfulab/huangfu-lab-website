# Flask Refactor Fix Report

**Date:** 2026-06-09 08:30:06
**Mode:** Auto
**Plan source:** `logs/flask_refactor_plan_2026-06-09_083006.md`
**Files modified:** `perturbseq_bp.py`
**Groups fixed:** 7 of 7 (3 bug-risk, 2 maintenance, 2 style/perf)

---

## Summary

All seven redundancy groups were resolved in a single pass against `perturbseq_bp.py`. The most impactful fixes eliminate a crash-on-sort SQL error (`n_ctg`), prevent silent divergence of the timepoint order list, and remove dead branching logic in the tooltip builder. Two maintenance items were extracted or aliased, and two style issues cleaned up module hygiene.

---

## Fixes Applied

### Group 1: `_SPARK_TP_ORDER` duplicates `_TP_ORDER` exactly — 🔴 Bug risk

**What was wrong:**
The timepoint order list `['ES_0h', 'DE_12h', 'DE_24h', 'DE_36h', 'DE_48h', 'DE_60h', 'DE_72h']` was literally copy-pasted to create `_SPARK_TP_ORDER` (line 380), a second constant with no drift guard. Any future addition to `_TP_ORDER` would leave sparkline rendering silently out of sync.

**What changed:**
```python
# Before (line 380)
_SPARK_TP_ORDER = ['ES_0h', 'DE_12h', 'DE_24h', 'DE_36h', 'DE_48h', 'DE_60h', 'DE_72h']

# After
_SPARK_TP_ORDER = _TP_ORDER
```

**Why this matters:**
Now sparklines are guaranteed to track `_TP_ORDER` automatically; no timepoint can be added to queries without also appearing in sparklines.

---

### Group 2: Ghost `n_ctg` sort column causes SQL error — 🔴 Bug risk

**What was wrong:**
`_query_tf_linked_genes_paged`'s `order_map` mapped column index 3 to `n_ctg`, a name that no longer exists in the SELECT (it was renamed to `n_linked_peaks`). Any frontend sort request on column 3 would raise a SQLite `no such column` error.

**What changed:**
```python
# Before
3: f'n_ctg {safe_dir}',

# After
3: f'n_linked_peaks {safe_dir}',
```

**Why this matters:**
Sorting by the linked-peaks column in the TF-linked-genes table now works instead of crashing.

---

### Group 3: Dead `if/else` with identical branches in `gene_page` — 🔴 Bug risk

**What was wrong:**
The `module_tooltips` builder in `gene_page` branched on `_r['source'] == 'hotspot_submodule'`, but both the `if` and `else` bodies were identical. The check added noise, misleadingly suggesting that submodule vs. supermodule tooltips are handled differently. `mt.source` was also being fetched and immediately discarded.

**What changed:**
```python
# Before
for _r in db.execute(
    f"SELECT mt.module_name, md.title, md.standard, mt.source ..."
).fetchall():
    if _r['source'] == 'hotspot_submodule':
        module_tooltips[_r['module_name']] = {'title': _r['title'], 'desc': _r['standard']}
    else:
        module_tooltips[_r['module_name']] = {'title': _r['title'], 'desc': _r['standard']}

# After
for _r in db.execute(
    f"SELECT mt.module_name, md.title, md.standard ..."
).fetchall():
    module_tooltips[_r['module_name']] = {'title': _r['title'], 'desc': _r['standard']}
```

**Why this matters:**
The dead branch no longer misleads future developers, and the query is one column lighter.

---

### Group 4: Duplicated `sources`/`datasets_str` parsing extracted to helper — 🟡 Maintenance burden

**What was wrong:**
An 8-line block parsing `sources_str` and `datasets_str` DB fields was copy-pasted identically into `_query_tf_linked_genes_paged` and `api_tf_linked_genes_all`. Any change to the `|||` encoding format required two edits.

**What changed:**
```python
# New helper added after rows_to_dicts
def _parse_sources_datasets(r) -> tuple:
    sources = sorted(set((r['sources_str'] or '').split(','))) if r['sources_str'] else []
    seen_ds: set = set()
    datasets = []
    for entry in (r['datasets_str'] or '').split(','):
        parts = entry.split('|||')
        if len(parts) == 4 and parts[0] not in seen_ds:
            seen_ds.add(parts[0])
            datasets.append({'id': parts[0], 'name': parts[1], 'cell_type': parts[2], 'source': parts[3]})
    return sources, datasets

# Both call sites replaced with:
sources, datasets = _parse_sources_datasets(r)
```

**Why this matters:**
The `|||`-encoding format is now maintained in one place; both API endpoints stay in sync automatically.

---

### Group 5: `BW_EXTRA_DIR` aliased to `BW_DIR` — 🟡 Maintenance burden

**What was wrong:**
`BW_EXTRA_DIR` was set to the same hardcoded path as `BW_DIR`. A future change to one would silently leave the other stale.

**What changed:**
```python
# Before
BW_EXTRA_DIR = str(Path(__file__).resolve().parent / "data" / "bw" / "atac")

# After
BW_EXTRA_DIR = BW_DIR
```

**Why this matters:**
`serve_bw_extra` will always serve from the same directory as `serve_bw`, tracked by a single source of truth.

---

### Group 6: Inline numpy/scipy imports moved to module level — 🔵 Style/perf

**What was wrong:**
`_coef_kde_data_db` imported `numpy` and `scipy.stats.gaussian_kde` inside the function body, running on every invocation and hiding dependencies from linters/IDEs.

**What changed:**
```python
# Before (inside _coef_kde_data_db)
import numpy as np
from scipy.stats import gaussian_kde

# After (top of file, with other imports)
import numpy as np
from scipy.stats import gaussian_kde
```

**Why this matters:**
Dependencies are visible at module level; type checkers and import-order linters now see them correctly.

---

### Group 7: Dead `get_tc_db` alias removed — 🔵 Style/perf

**What was wrong:**
`get_tc_db` was a one-line wrapper around `get_db` with no callers inside the file and no external imports.

**What changed:**
```python
# Removed entirely:
def get_tc_db():
    return get_db()
```

**Why this matters:**
Dead code is gone; developers no longer wonder whether `get_tc_db` and `get_db` behave differently.

---

## Groups Not Fixed

None — all groups from the plan were resolved.

---

## File Change Index

| File | Nature of change |
|------|-----------------|
| `perturbseq_bp.py` | Added `import numpy as np` and `from scipy.stats import gaussian_kde` at top |
| `perturbseq_bp.py` | `BW_EXTRA_DIR = BW_DIR` (was duplicate path literal) |
| `perturbseq_bp.py` | Removed dead `get_tc_db()` alias |
| `perturbseq_bp.py` | Added `_parse_sources_datasets(r)` helper after `rows_to_dicts` |
| `perturbseq_bp.py` | `_SPARK_TP_ORDER = _TP_ORDER` (was duplicate list literal) |
| `perturbseq_bp.py` | `order_map[3]`: `n_ctg` → `n_linked_peaks` in `_query_tf_linked_genes_paged` |
| `perturbseq_bp.py` | Collapsed dead `if/else` in `gene_page` module_tooltips builder; removed `mt.source` from SELECT |
| `perturbseq_bp.py` | Both `data` loops in `_query_tf_linked_genes_paged` and `api_tf_linked_genes_all` use `_parse_sources_datasets(r)` |
| `perturbseq_bp.py` | Removed inline imports from `_coef_kde_data_db` |
