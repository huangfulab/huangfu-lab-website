# Refactor Notes: Code Duplication & Consolidation Opportunities

> Audit date: 2026-06-03. No code was changed; this document describes what *could* be improved.

---

## 1. Python: Repeated SQL Queries

### 1.1 Navigation lists (TFs + supermodules) — HIGH

`_nav_lists()` (line ~754) already encapsulates the two queries that load the sidebar navigation lists:

```python
tfs = [r[0] for r in db.execute(
    "SELECT gene_name FROM gene_table WHERE in_perturbation_library=1 ORDER BY gene_name")]
modules = [r[0] for r in db.execute(
    "SELECT module_name FROM module_table "
    "WHERE source='hotspot_supermodule' AND module_name != 'unassigned' ORDER BY module_name")]
```

**Where duplicated:** Inline in the index route (`perturbseq_bp.py` ~line 437–441) and again inside `_nav_lists()` (~line 754–758). Every route that renders a full page calls `_nav_lists()`, but the index route re-queries manually.

**Risk:** If the query is ever changed (e.g. adding a filter or join), the index route will silently diverge.

**Fix:** Replace the inline queries in the index route with a call to `_nav_lists(db)`.

---

### 1.2 Gene module membership — HIGH

The submodule/supermodule membership of a gene is resolved with:

```python
mem_rows = db.execute("""
    SELECT m.module_name, m.source
    FROM gene_module_table gm JOIN module_table m ON gm.module_id = m.module_id
    WHERE gm.gene_id=?
""", (gene_id,)).fetchall()
submodule = next((r['module_name'] for r in mem_rows if r['source'] == 'hotspot_submodule'), None)
supermodule = next((r['module_name'] for r in mem_rows if r['source'] == 'hotspot_supermodule'), None)
```

**Where duplicated:** `gene_page` (~line 1232–1244), `api_gene` (~line 1872–1879), `api_gene_tooltip` (~line 2123–2133). The `api_gene` version skips the `!= 'unassigned'` guard that `gene_page` applies.

**Risk:** The `unassigned` guard is inconsistently applied — some callers filter it out, others expose it to the frontend.

**Fix:** Extract to `_get_gene_membership(db, gene_id) -> dict` that always normalises `'unassigned'` to `None`.

---

### 1.3 Module genes list — HIGH

Fetching the gene members of a module:

```python
genes = [r[0] for r in db.execute(
    "SELECT gt.gene_name FROM gene_module_table gm "
    "JOIN gene_table gt ON gm.gene_id = gt.gene_id WHERE gm.module_id=?",
    (module_id,))]
```

**Where duplicated:** `_super_page` (~line 774–777), `_sub_page` (~line 892–896), `_gc_page` (~line 1012–1017), `api_module_detail` (~line 2033–2035), `api_module_genes` (~line 2073–2075). Some versions include `ORDER BY gt.gene_name`, others don't.

**Risk:** The inconsistent ordering means the gene list can appear in different orders on different pages for the same module.

**Fix:** Extract to `_get_module_genes(db, module_id, ordered=True) -> list[str]` and use everywhere.

---

### 1.4 Module description — MEDIUM

Fetching the `title` and `standard` description for a module:

```python
mod_desc_row = db.execute(
    "SELECT title, standard FROM module_description WHERE module_id=?",
    (module_id,)).fetchone()
```

**Where duplicated:** `_super_page` (~line 829–830), `_sub_page` (~line 937–938), `_gc_page` (~line 1035–1037), `api_module_tooltip` (~line 2088–2090).

**Risk:** Low — the query is simple and stable — but it still adds noise.

**Fix:** Extract to `_get_module_description(db, module_id) -> dict | None`.

---

### 1.5 GO enrichment highlight terms — MEDIUM

Top GO:BP terms for a module card header:

```python
highlight_terms = [r[0] for r in db.execute(
    "SELECT term_name FROM go_module_enrichment "
    "WHERE module_id=? AND source='GO:BP' ORDER BY p_value LIMIT 5",
    (module_id,)).fetchall()]
```

**Where duplicated:** `_super_page` (~line 831–834, LIMIT 5), `_sub_page` (~line 914–917, LIMIT 4), `_gc_page` (~line 1039–1042, LIMIT 5).

**Risk:** The limit differs (4 vs 5) without explanation. A future change to the display logic may only update one call site.

**Fix:** Extract to `_get_module_go_highlights(db, module_id, limit=5) -> list[str]` and pass an explicit limit where needed.

---

### 1.6 TF perturbation + binding merge pattern — MEDIUM

Building the merged TF regulation table (perturbation NES + binding odds ratio) appears six times across route handlers and API endpoints with only minor differences in the `module_collection` / `gene_set_collection` filter values:

```python
pert_rows = rows_to_dicts(db.execute(
    "SELECT module, mean_NES, n_grnas AS n_sig_gRNA, min_padj AS padj "
    "FROM gsea_tf_table WHERE gene_name=? AND module_collection=? "
    "AND gene_set_collection=?", (...)))
bind_by_mod = _bind_edges_for_tf(db, gene_name, collection)
tf_reg = _merge_pert_bind_edges(pert_rows, bind_by_mod, id_key='module')
```

**Where duplicated:** `_super_page` (~805–821), `_sub_page` (~918–923), `_gc_page` (~972–1010), `gene_page` (~1286–1299 for supermodule, ~1360–1395 for submodule), `tf_module_link_page` (~1752–1757).

**Risk:** If the `gsea_tf_table` schema changes (e.g. column rename), every call site needs updating.

**Fix:** Extract to `_get_tf_regulation(db, gene_name, module_collection, gene_set_collection) -> list[dict]` that calls `_bind_edges_for_tf` and `_merge_pert_bind_edges` internally.

---

### 1.7 Gene TF status classification — HIGH

The three-way classification (Active TF / TF / Gene) is performed inline in three route handlers despite a `_build_gene_data_list()` helper already existing:

```python
for g in genes:
    if g in perturbed:
        tf_status = 'Active TF'
    elif g in lambert:
        tf_status = 'TF'
    else:
        tf_status = 'Gene'
```

**Where duplicated (inline):** `_sub_page` (~900–908), `_gc_page` (~1022–1029), `api_module_detail` (~2040–2047).

**Where the helper exists but is not used:** `_build_gene_data_list()` at ~line 714 is called from `go_page` and a couple of other places, but not from the three above.

**Risk:** If the classification logic changes (e.g. a fourth status category), the three inline copies will silently produce different output from pages that use the helper.

**Fix:** Replace all three inline blocks with calls to `_build_gene_data_list()`, adjusting the `name_key` / `pub_key` arguments as needed.

---

## 2. Python: Long Route Handlers

### 2.1 `gene_page` (~262 lines, ~line 1212–1474)

This single route handler performs:
- Gene ID lookup + 404 guard
- Module membership resolution (submodule + supermodule)
- TF vs non-TF branching (15+ DB queries for TF path)
- Expression data (timecourse + bulk)
- Regulatory element genomic data
- Tooltip and chart payload assembly

**Risk:** The function is difficult to unit-test and easy to introduce regressions in when adding features to one branch (TF vs non-TF).

**Fix:** Split into private helpers: `_get_gene_basic_info`, `_get_gene_modules`, `_get_gene_tf_data`, `_get_gene_expression`, then call them from a thin route handler.

---

### 2.2 `_query_element_genes` (~78 lines, ~line 2449–2527) and `_query_tf_gene_link` (~82 lines, ~line 2305–2387)

Both functions combine complex CTE queries with manual deduplication loops. The deduplication logic (building `by_peak` dicts with `_seen` bookkeeping keys) in `_query_tf_gene_link` is particularly hard to follow.

**Fix:** Split each into: (a) raw SQL fetch, (b) deduplicate/group, (c) format output.

---

## 3. Templates: Repeated Structure

### 3.1 Column definitions for data tables — LOW

The column metadata dict (passed to `makeFilterSortTable`) is copy-pasted with minor differences across three templates:

| Template | Variable | Key difference |
|---|---|---|
| `tf.html` (~line 3–11) | `cols_tf` | 6 columns |
| `go_term.html` (~line 5–13) | `cols_tf` | 7 columns (adds `odds_ratio`) |
| `module.html` (~line 128–136) | `cols_full` | 7 columns (same as go_term) |

**Risk:** If a shared column (e.g. `mean_NES` label or format) changes, each template must be updated separately.

**Fix:** Define canonical column lists in Python as module-level constants and pass them in template context, or create a shared `_table_columns.html` include.

---

### 3.2 Filter bar: mixed macro and inline HTML — MEDIUM

Three approaches currently exist for the evidence/direction filter bar:

1. **Inline HTML** (`tf.html` ~line 125–134): hardcoded `<button>` tags
2. **Macro call** (`go_term.html` ~line 89): `{{ filter_bar(...) }}` from `_network_macros.html`
3. **Split inline HTML** (`module.html`): Evidence and direction bars coded separately

**Risk:** CSS class names for the buttons differ slightly between implementations, making JS event listeners fragile if one template is updated.

**Fix:** Standardise on the `filter_bar` macro from `_network_macros.html` and extend it to support the split evidence/direction variant.

---

### 3.3 `landing_split` macro not used in `gene.html` — LOW

The `_network_macros.html` macro `landing_split` (Cytoscape panel + side table layout) is imported and used in `tf.html`, `go_term.html`, and `module.html`, but `gene.html` builds the same two-column layout with raw HTML.

**Risk:** Visual inconsistencies if the layout is updated in the macro but not in `gene.html`.

**Fix:** Import `_network_macros.html` in `gene.html` and replace the raw layout with `landing_split`.

---

### 3.4 Hero / page header structure — LOW

Each landing page defines a page hero with a name and badge:

```html
<!-- tf.html -->
<div class="tf-hero"> <h1 class="tf-symbol">{{ tf }}</h1> <span class="tf-badge">Active TF</span> </div>

<!-- go_term.html -->
<div class="go-hero"> <h1 class="go-name">{{ go_name }}</h1> <span class="go-badge">GO Term</span> </div>

<!-- module.html -->
<div class="mod-hero"> <span class="mod-symbol">{{ mod }}</span> <span class="mod-badge">Supermodule</span> </div>
```

The HTML structure and intent are identical; only CSS class prefixes and label text differ.

**Risk:** Low currently, but divergence will grow as pages are independently restyled.

**Fix:** Create a shared `_hero_header.html` macro: `{% macro hero(name, badge, color_class='') %}`.

---

### 3.5 Gene + module tooltip system — MEDIUM

`_tooltip.html` defines two nearly-identical tooltip components (lines ~2–107 for genes, ~109–234 for modules). Each implements:

1. `mouseover` detection via `closest('[href]')`
2. `fetch` to the relevant API endpoint
3. In-memory cache (`Map`)
4. HTML render into a shared `<div class="gene-tooltip">`

The two implementations share ~80% of their code structure.

**Risk:** Bug fixes or styling changes to the tooltip mechanism must be applied twice.

**Fix:** Unify into a single configurable tooltip component, parameterised by `linkSelector`, `apiEndpoint`, and `renderFn`.

---

## 4. JavaScript: Repeated Patterns

### 4.1 Two similar table utilities — LOW

`tf-network.js` contains two table-rendering functions with overlapping logic:

| Function | Lines | Capabilities |
|---|---|---|
| `makeFilterSortTable` | ~153–273 | Filter by evidence/direction, sort, search, page |
| `makeEnrichTable` | ~308–378 | Sort, page (no evidence/direction filter) |

Both share: column rendering, sort toggling, and page-navigation logic.

**Risk:** Pagination or sort logic fixed in one will silently remain broken in the other.

**Fix:** Extract shared primitives (renderRows, sortToggle, paginate) into private helpers; have both functions call them.

---

### 4.2 Cytoscape initialisation repeated per page — MEDIUM

`initNetwork()` in `tf-network.js` provides a base, but each page's inline `<script>` block re-implements: layout configuration, `tap` event handlers, and `mouseover` tooltip wiring. The three pages (`tf.html`, `go_term.html`, `module.html`) each have 40–80 lines of near-identical Cytoscape setup.

**Risk:** A fix to the tap-handler logic (e.g. stopping propagation) must be applied in three places.

**Fix:** Extend `TFNetwork` with `initFullNetwork(container, nodes, edges, opts)` that accepts callback hooks for tap/hover behaviour rather than hardcoding them per page.

---

### 4.3 `MODULE_COLORS` defined in four places — MEDIUM

The supermodule colour palette is hardcoded in:

1. `perturbseq_bp.py` (Python constant, injected into templates)
2. `perturbseq/index.html` (inline JS object)
3. `static/js/tf-network.js` (JS constant)
4. `perturbseq/gene.html` (inline JS object)

(This is already noted in `CODE_AUDIT.md` but included here for completeness.)

**Risk:** Adding or recolouring a supermodule requires four synchronised edits.

**Fix:** Store the palette in a single JSON file (or as a DB table); load it once in Python and inject into all templates via `g` or a context processor; have `tf-network.js` consume it from the injected `window.MODULE_COLORS`.

---

### 4.4 TF status badge rendered in JavaScript inline strings — LOW

Several templates generate TF-status badge HTML by string concatenation inside JavaScript:

```javascript
const badge = r.tf_status === 'Active TF'
  ? `<span class="intro-tag gene-tag" style="background:var(--accent)">Active TF</span>`
  : r.tf_status === 'TF'
  ? `<span class="intro-tag gene-tag" style="background:#888">TF</span>`
  : '';
```

This pattern appears in `go_term.html`, `module.html`, and `gene.html`.

**Risk:** CSS class names and inline styles in JS strings are invisible to linters; a class rename will silently break badges.

**Fix:** Move badge generation into a shared JS function `tfStatusBadge(status)` inside `tf-network.js` and call it from all table renderers.

---

## 5. Priority Summary

| # | Issue | File(s) | Impact | Effort |
|---|---|---|---|---|
| 1 | Gene TF status inline vs helper | `perturbseq_bp.py` | High | Low |
| 2 | Gene module membership helper | `perturbseq_bp.py` | High | Low |
| 3 | Module genes query | `perturbseq_bp.py` | High | Low |
| 4 | Index route duplicates `_nav_lists` | `perturbseq_bp.py` | High | Trivial |
| 5 | `MODULE_COLORS` in 4 places | py + 3 templates | High | Medium |
| 6 | TF perturbation+binding merge | `perturbseq_bp.py` | Medium | Low |
| 7 | Module description helper | `perturbseq_bp.py` | Medium | Low |
| 8 | GO highlight terms helper | `perturbseq_bp.py` | Medium | Low |
| 9 | Gene+module tooltip unification | `_tooltip.html` | Medium | Medium |
| 10 | Filter bar macro vs inline HTML | `tf.html`, `go_term.html`, `module.html` | Medium | Low |
| 11 | Cytoscape init per-page repetition | `tf.html`, `go_term.html`, `module.html` | Medium | High |
| 12 | `gene_page` 262-line handler | `perturbseq_bp.py` | Medium | Medium |
| 13 | Two JS table utilities overlap | `tf-network.js` | Low | Medium |
| 14 | `landing_split` not used in gene page | `gene.html` | Low | Low |
| 15 | Column definitions copy-pasted | `tf.html`, `go_term.html`, `module.html` | Low | Low |
| 16 | TF status badge in JS strings | `go_term.html`, `module.html`, `gene.html` | Low | Low |
| 17 | Hero header structure per page | all landing pages | Low | Low |
