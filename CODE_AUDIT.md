# Webapp Code Audit

> Status: findings only — no edits made. Ordered by severity/impact.

---

## 1. MODULE COLORS — quadruple duplication (CRITICAL)

`MODULE_COLORS` is hardcoded in four separate places that must be kept in sync:

| Location | Form |
|---|---|
| `perturbseq_bp.py:23–28` | Python dict |
| `templates/.../index.html:133–138` | JS literal |
| `static/js/tf-network.js:9–14` | JS literal |
| `templates/.../gene.html:315` | single color derived from it |

**Fix**: Keep only the Python dict; inject into all templates as
`const MODULE_COLORS = {{ module_colors \| tojson }};`.

---

## 2. `drawExprChart()` — duplicated D3 function

Identical D3 v7 line chart (margin/scale/axis/area/tooltip) is copy-pasted in:

- `templates/perturbseq/tf_gene_link.html:406–475`
- `templates/perturbseq/gene.html:320–421` (same pattern, slightly different wrapper)

Both use `scalePoint`, `scaleLinear`, `curveCatmullRom`, the fixed-div tooltip, and circle
hover handlers. Only the colour argument differs.

**Fix**: Extract to `static/js/expr-chart.js`, load once in `base.html` (or `extra_style`),
call from both templates.

---

## 3. Chart axis CSS — conflicting duplicate rules

`.axis text` font-size differs between the global stylesheet and an inline override:

- `static/css/perturbseq.css:~864` → `font-size: 10px`
- `templates/perturbseq/tf_gene_link.html:100` → `font-size: 9px` (inline `<style>`)

The following rules are also duplicated between the CSS file and `tf_gene_link.html`'s inline block:
`.axis path`, `.axis line`, `.grid line`, `.grid path`, `.area-fill`, `.mean-line`.

**Fix**: Remove the inline chart-axis rules from `tf_gene_link.html`; reconcile the
font-size to a single value in `perturbseq.css`.

---

## 4. Database connection inconsistency

Three different patterns are used across `perturbseq_bp.py`:

| Pattern | Lines | Problem |
|---|---|---|
| `get_db()` — Flask `g`, `sqlite3.Row` factory | 210, 330, 377, … | Correct pattern |
| `sqlite3.connect()` + manual `.close()` | 185–190, 198–203 | No `row_factory`; rows are plain tuples accessed by index; no context manager (connection leaks on exception) |
| `get_tc_db()` named inconsistently | called `tc` (l.465, l.555) vs `tc_db` (l.706) | Style inconsistency |

**Fix**: Replace bare `sqlite3.connect()` calls with `get_db()` (or a context manager);
standardise timecourse variable name to `tc`.

---

## 5. Inconsistent row access (tuple index vs dict key)

`get_db()` sets `sqlite3.Row` so columns should be accessed by name, but several
call sites still use positional indexing:

- `perturbseq_bp.py:241–242` — `[r[0] for r in db.execute(...)]`
- Compare `perturbseq_bp.py:362` — `dict(node_row)` (correct)

**Fix**: Standardise to `r["column_name"]` everywhere.

---

## 6. `json.loads()` without error handling

Two places deserialise DB-stored JSON without a try/except:

- `perturbseq_bp.py:337` — `json.loads(desc["highlight_terms"])`
- `perturbseq_bp.py:517–518` — similar pattern

A malformed stored value raises an unhandled exception → HTTP 500.

**Fix**: Wrap in `try/except json.JSONDecodeError` and return a sensible default.

---

## 7. `link_type_meta` Jinja dict — defined per-template

`link_type_meta` is a local `{% set %}` in `tf_gene_link.html:198–203`.  Any future
template that needs the same colour/label mapping will re-define it.

**Fix**: Move into a Jinja macro file (e.g. `templates/perturbseq/macros.html`) or pass
from the route as a template variable; share across templates via `{% from … import %}`.

---

## 8. Hardcoded URLs — should use `url_for()`

More than a dozen places build URLs with string concatenation instead of `url_for()`.
If the blueprint prefix or route ever changes, these all silently break:

- `templates/perturbseq/gene.html:126, 131, 136, 137` — `/perturbseq/module/…`
- `templates/perturbseq/base.html:88–89` — `/perturbseq/module/`, `/perturbseq/gene/`
- `templates/perturbseq/tf_gene_link.html:219, 226` — `/perturbseq/gene/{{ tf }}`, `/perturbseq/gene/{{ gene }}`
- `perturbseq_bp.py:201, 252, 631, 635` (Python string builds)

**Fix**: Replace with `url_for('perturbseq.gene_page', gene=tf)` etc.

---

## 9. Cytoscape instances leaked to `window`

Three separate Cytoscape network graphs in `gene.html` store their instances globally:

- `gene.html:463` — `window._tfModCy = cy;`
- `gene.html:547` — `window._tfDevCy = cy;`
- `gene.html:708` — `window._tfSubCy = cy;`

There is no teardown / `cy.destroy()` on navigation, which can leak memory in
single-page-like navigation patterns.

Also `window.clearFilters` in `tf_gene_link.html:661` pollutes the global scope
unnecessarily (it is only used from an `onclick` attribute on one button).

**Fix**: Scope Cytoscape instances to their initialising IIFE and only expose the minimum
needed. Replace `onclick="clearFilters()"` with a proper `addEventListener`.

---

## 10. `_build_gene_data_list()` — not used consistently

The helper at `perturbseq_bp.py:307–318` exists precisely to build gene data objects
with TF-status flags, but the same logic is re-inlined at `perturbseq_bp.py:387–394`.

**Fix**: Replace the inlined copy with a call to the existing helper.

---

## 11. D3 tooltip div — created per chart, never removed

`drawExprChart()` does:
```js
const tip = d3.select('body').append('div')…
```
Each call appends a new `<div>` to `<body>`. If `drawExprChart` is called more than once
(e.g. on resize/redraw), tooltip divs accumulate.

**Fix**: Create the tooltip div once (shared across charts) or select-or-create with
`d3.select('#chart-tooltip')`.

---

## 12. D3 / IGV loaded in `extra_style` block

`tf_gene_link.html:5–7` loads `<script>` and `<link>` tags inside `{% block extra_style %}`.
Semantically this block is for CSS; scripts loaded here appear in `<head>` but before the
closing `</head>` which is fine, however it conflates two concerns.

**Fix**: Add an `{% block extra_head %}` block to `base.html` for arbitrary head content;
keep `extra_style` for `<link rel="stylesheet">` only.

---

## 13. Potentially dead CSS classes

The following classes are defined in `perturbseq.css` but appear to have no matching HTML
in any current template (requires final grep to confirm):

- `.enrich-bar-cell`, `.enrich-bar`, `.enrich-source-filters`, `.enrich-src-btn` (~lines 488–540)
- `.tf-type-active`, `.tf-type-tf`

**Fix**: Confirm with `grep -r "enrich-bar\|tf-type-" templates/`; remove if truly unused.

---

## 14. `_TF_GENE_LINK_DUMMY` — route always returns dummy data

`perturbseq_bp.py:1227–1233`: The `/link/<tf>/<gene>` route ignores the `tf` and `gene`
URL parameters entirely and always renders the hard-coded dummy dict. No 404 is returned
for non-existent pairs.

This is intentional for now but should be tracked as a known gap until the DB tables
(`tf_gene_links`, `tf_element_motif`, etc.) are populated.

---

## Summary table

| # | Category | Files affected | Severity |
|---|---|---|---|
| 1 | MODULE_COLORS 4× duplication | bp.py, index.html, tf-network.js, gene.html | High |
| 2 | `drawExprChart` copy-pasted | gene.html, tf_gene_link.html | High |
| 3 | Chart CSS duplicated / conflicting | perturbseq.css, tf_gene_link.html | Medium |
| 4 | DB connection pattern mismatch | perturbseq_bp.py | Medium |
| 5 | Row access: index vs name | perturbseq_bp.py | Medium |
| 6 | `json.loads` without try/except | perturbseq_bp.py | Medium |
| 7 | `link_type_meta` per-template | tf_gene_link.html | Low-Medium |
| 8 | Hardcoded URLs (12+ sites) | gene.html, base.html, tf_gene_link.html, bp.py | Low-Medium |
| 9 | Cytoscape / clearFilters on `window` | gene.html, tf_gene_link.html | Low |
| 10 | `_build_gene_data_list` not used consistently | perturbseq_bp.py | Low |
| 11 | D3 tooltip div accumulates | gene.html, tf_gene_link.html | Low |
| 12 | Scripts in `extra_style` block | tf_gene_link.html, base.html | Low |
| 13 | Possibly dead CSS | perturbseq.css | Low |
| 14 | `/link` route ignores URL params | perturbseq_bp.py | Known gap |
