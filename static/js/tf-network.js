/**
 * Shared utilities for TF regulatory network displays.
 * Used across module, TF, and gene pages.
 */
window.TFNetwork = (function () {
  'use strict';

  // ── Module color palette — injected by server via base.html; fallback for safety ──
  const MODULE_COLORS = window.MODULE_COLORS || {
    'DE-1':'#4E79A7','DE-2':'#A0CBE8','DE-3':'#F28E2B',
    'DE-4':'#FFBE7D','DE-5':'#59A14F','DE-6':'#8CD17D','DE-7':'#B6992D',
    'DE-8':'#F1CE63','DE-9':'#499894','DE-10':'#86BCB6','DE-11':'#E15759',
    'DE-12':'#FF9D9A',
  };

  // ── Cytoscape style presets ───────────────────────────────────────────────

  /** Full edge styles: perturbation+binding, perturbation-only, binding-only */
  const EDGE_STYLES_FULL = [
    { selector:'edge[direction="up"][evidence="both"]',           style:{'line-color':'#ec4899','width':2.5,'curve-style':'bezier','target-arrow-color':'#ec4899','target-arrow-shape':'triangle'} },
    { selector:'edge[direction="down"][evidence="both"]',         style:{'line-color':'#3b82f6','width':2.5,'curve-style':'bezier','target-arrow-color':'#3b82f6','target-arrow-shape':'triangle'} },
    { selector:'edge[direction="up"][evidence="perturbation"]',   style:{'line-color':'#ec4899','width':2,'line-style':'dashed','curve-style':'bezier','target-arrow-color':'#ec4899','target-arrow-shape':'triangle'} },
    { selector:'edge[direction="down"][evidence="perturbation"]', style:{'line-color':'#3b82f6','width':2,'line-style':'dashed','curve-style':'bezier','target-arrow-color':'#3b82f6','target-arrow-shape':'triangle'} },
    { selector:'edge[direction="bind"]',                          style:{'line-color':'#9098b8','width':1.5,'line-style':'dashed','curve-style':'bezier','target-arrow-color':'#9098b8','target-arrow-shape':'triangle'} },
    { selector:'edge[evidence="binding"]',                        style:{'line-color':'#9098b8','width':1.5,'line-style':'dashed','curve-style':'bezier','target-arrow-color':'#9098b8','target-arrow-shape':'triangle'} },
  ];

  /** GC/PC edge styles: only up/down (no evidence distinction) */
  const EDGE_STYLES_GC = [
    { selector:'edge[direction="up"]',   style:{'line-color':'#ec4899','width':2,'line-style':'dashed','curve-style':'bezier','target-arrow-color':'#ec4899','target-arrow-shape':'triangle'} },
    { selector:'edge[direction="down"]', style:{'line-color':'#3b82f6','width':2,'line-style':'dashed','curve-style':'bezier','target-arrow-color':'#3b82f6','target-arrow-shape':'triangle'} },
  ];

  /** Common node style for a central "hub" node (module or TF) */
  function hubNodeStyle(opts) {
    return {
      'background-color': opts.colorAttr || '#4a56b0',
      'label': 'data(label)',
      'font-size': opts.fontSize || '12px',
      'color': '#1e2448',
      'text-valign': 'bottom',
      'text-margin-y': 5,
      'width': opts.size || 50,
      'height': opts.size || 50,
      'shape': opts.shape || 'round-rectangle',
      'border-width': 1,
      'border-color': 'rgba(0,0,0,0.15)',
    };
  }

  /** Common node style for peripheral nodes (TFs, modules, submodules) */
  function peripheralNodeStyle(opts) {
    return {
      'background-color': opts.colorAttr || '#4a56b0',
      'label': 'data(label)',
      'font-size': opts.fontSize || '9px',
      'color': opts.textColor || '#1e2448',
      'text-valign': opts.textValign || 'bottom',
      'text-halign': opts.textHalign || 'center',
      'text-margin-y': opts.textMarginY != null ? opts.textMarginY : 4,
      'width': opts.size || 22,
      'height': opts.size || 22,
      'shape': opts.shape || 'ellipse',
      'border-width': opts.borderWidth != null ? opts.borderWidth : 0,
      'border-color': opts.borderColor || 'rgba(0,0,0,0.2)',
    };
  }

  // ── Position utility ──────────────────────────────────────────────────────

  /**
   * Arrange items in a sector (arc segment) around the origin.
   * @param {string[]} ids        - Node IDs to position
   * @param {number}   centerDeg  - Center angle in degrees (0=right, 90=down, 180=left)
   * @param {number}   spreadDeg  - Angular spread in degrees
   * @param {number}   rStart     - Starting radius
   * @param {number}   rStep      - Radius increment per row
   * @param {number}   perRow     - Max items per row before wrapping
   * @returns {Object} Map of id → {x, y}
   */
  function placeSector(ids, centerDeg, spreadDeg, rStart, rStep, perRow) {
    const pos = {};
    ids.forEach((id, i) => {
      const row = Math.floor(i / perRow);
      const col = i % perRow;
      const r = rStart + row * rStep;
      const count = Math.min(perRow, ids.length - row * perRow);
      const half = count > 1 ? spreadDeg / 2 : 0;
      const step = count > 1 ? spreadDeg / (count - 1) : 0;
      const angle = (centerDeg - half + col * step) * Math.PI / 180;
      pos[id] = { x: r * Math.cos(angle), y: r * Math.sin(angle) };
    });
    return pos;
  }

  // ── Cytoscape initializer ─────────────────────────────────────────────────

  /**
   * Initialize a Cytoscape network instance.
   * @param {Object} opts
   * @param {string|HTMLElement} opts.container - ID or element
   * @param {Array}  opts.elements  - Cytoscape elements array
   * @param {Array}  [opts.style]   - Full cytoscape style array (overrides preset)
   * @param {string} [opts.preset]  - 'full' | 'gc' (used if opts.style not provided)
   * @param {Array}  [opts.nodeStyles] - Additional node style selectors
   * @param {Object} [opts.layout]  - Cytoscape layout config (default: preset)
   * @param {Function} [opts.onTapNode] - Callback(nodeId, nodeData) on node tap
   * @param {string} [opts.tapSelector] - Selector for tap events (default: 'node')
   * @returns {Object} Cytoscape instance
   */
  function initNetwork(opts) {
    const container = typeof opts.container === 'string'
      ? document.getElementById(opts.container)
      : opts.container;

    const edgeStyles = opts.style ? [] : (opts.preset === 'gc' ? EDGE_STYLES_GC : EDGE_STYLES_FULL);
    const style = opts.style || [...(opts.nodeStyles || []), ...edgeStyles];

    const cy = cytoscape({
      container: container,
      elements: opts.elements,
      minZoom: opts.minZoom != null ? opts.minZoom : 0.15,
      wheelSensitivity: opts.wheelSensitivity || 0.25,
      style: style,
      layout: opts.layout || { name: 'preset', animate: false },
    });

    if (opts.onTapNode) {
      const sel = opts.tapSelector || 'node';
      cy.on('tap', sel, evt => opts.onTapNode(evt.target.id(), evt.target.data()));
    }

    return cy;
  }

  // ── Filter + Sort Table ───────────────────────────────────────────────────

  /**
   * Generic filter + sort table with optional filter bar and search.
   *
   * @param {HTMLElement} table     - The <table> element
   * @param {Array}       allData   - Full data array
   * @param {Function}    renderRow - (rowData) => HTMLTableRowElement
   * @param {Object}      [opts]    - Options
   * @param {string}      [opts.barId]     - ID of filter bar element (buttons with data-filter or data-ev/data-dir)
   * @param {string}      [opts.countId]   - ID of count display element
   * @param {string}      [opts.searchId]  - ID of search <input> element
   * @param {string}      [opts.searchKey] - Data key to search on (default: auto-detect 'tf' or 'module')
   * @param {Function}    [opts.onChange]  - Callback after each refresh
   * @param {Function}    [opts.filterFn]  - Custom filter function (data, {evFilter, dirFilter, searchQ}) => filteredData
   * @returns {Object}    Control object with { refresh(), getFiltered() }
   */
  function makeFilterSortTable(table, allData, renderRow, opts) {
    opts = opts || {};
    const headers = table.querySelectorAll('th.sortable');
    let sortKey = null, sortDir = 'asc';
    let evFilter = 'all', dirFilter = 'all', searchQ = '';
    const PAGE_SIZE = opts.pageSize || 10;
    let page = 0;

    // Inject a pagination bar immediately after the table
    const pagEl = document.createElement('div');
    pagEl.className = 'table-pag';
    pagEl.style.display = 'none';
    table.parentNode.insertBefore(pagEl, table.nextSibling);

    // Auto-detect search key
    const searchKey = opts.searchKey || (allData[0] && 'tf' in allData[0] ? 'tf' : 'module');

    function applyFilter(d) {
      if (opts.filterFn) return opts.filterFn(d, { evFilter, dirFilter, searchQ });
      let r = d;
      if (evFilter === 'perturbation') r = r.filter(x => x.evidence === 'perturbation' || x.evidence === 'both');
      else if (evFilter === 'both')    r = r.filter(x => x.evidence === 'both');
      else if (evFilter === 'binding') r = r.filter(x => x.evidence === 'binding' || x.evidence === 'both');
      if (dirFilter === 'up')          r = r.filter(x => x.direction === 'Up');
      else if (dirFilter === 'down')   r = r.filter(x => x.direction === 'Down');
      if (searchQ) r = r.filter(x => {
        const v = x[searchKey] || '';
        return v.toLowerCase().includes(searchQ);
      });
      return r;
    }

    function sortData(d) {
      if (!sortKey) return d;
      const th = [...headers].find(h => h.dataset.key === sortKey);
      const type = th ? th.dataset.type : 'string';
      return [...d].sort((a, b) => {
        let va = a[sortKey], vb = b[sortKey];
        if (va == null) return 1;
        if (vb == null) return -1;
        if (type === 'number') return sortDir === 'asc' ? va - vb : vb - va;
        return sortDir === 'asc' ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
      });
    }

    function refresh() {
      const f = applyFilter(allData);
      const s = sortData(f);
      const tb = table.querySelector('tbody');
      tb.innerHTML = '';
      const totalPages = Math.max(1, Math.ceil(s.length / PAGE_SIZE));
      page = Math.min(page, totalPages - 1);
      s.slice(page * PAGE_SIZE, (page + 1) * PAGE_SIZE).forEach(r => tb.appendChild(renderRow(r)));
      pagEl.style.display = totalPages <= 1 ? 'none' : 'flex';
      pagEl.innerHTML = '<button class="table-pag-btn"' + (page === 0 ? ' disabled' : '') + '>← Prev</button>'
        + '<span>' + (page + 1) + ' / ' + totalPages + '</span>'
        + '<button class="table-pag-btn"' + (page >= totalPages - 1 ? ' disabled' : '') + '>Next →</button>';
      pagEl.querySelectorAll('.table-pag-btn').forEach((btn, i) => {
        btn.addEventListener('click', () => { page += (i === 0 ? -1 : 1); refresh(); });
      });
      if (opts.countId) {
        const el = document.getElementById(opts.countId);
        if (el) el.textContent = `${s.length} of ${allData.length}`;
      }
      if (opts.onChange) opts.onChange();
    }

    // Sortable headers
    headers.forEach(th => th.addEventListener('click', () => {
      const key = th.dataset.key, type = th.dataset.type;
      if (sortKey === key) sortDir = sortDir === 'asc' ? 'desc' : 'asc';
      else { sortKey = key; sortDir = type === 'number' ? 'desc' : 'asc'; }
      headers.forEach(h => h.classList.remove('sort-asc', 'sort-desc'));
      th.classList.add(sortDir === 'asc' ? 'sort-asc' : 'sort-desc');
      page = 0; refresh();
    }));

    // Filter bar — supports two modes:
    // 1) Simple: buttons with data-filter (tf.html style)
    // 2) Split: buttons with data-ev + data-dir (module.html style)
    if (opts.barId) {
      const bar = document.getElementById(opts.barId);
      if (bar) {
        const simpleBtns = bar.querySelectorAll('button[data-filter]');
        const evBtns = bar.querySelectorAll('button[data-ev]');
        const dirBtns = bar.querySelectorAll('button[data-dir]');

        if (evBtns.length) {
          // Split mode (module pages)
          function syncDirDisabled() {
            const bindingOnly = evFilter === 'binding';
            dirBtns.forEach(b => { b.disabled = bindingOnly; b.classList.toggle('filter-disabled', bindingOnly); });
            const dirLabel = bar.querySelector('.filter-label[id$="-dir-label"]') || document.getElementById('dir-label');
            if (dirLabel) dirLabel.classList.toggle('filter-disabled', bindingOnly);
            if (bindingOnly) { dirFilter = 'all'; dirBtns.forEach(b => b.classList.toggle('active', b.dataset.dir === 'all')); }
          }
          evBtns.forEach(btn => btn.addEventListener('click', () => {
            evBtns.forEach(b => b.classList.remove('active'));
            btn.classList.add('active'); evFilter = btn.dataset.ev; syncDirDisabled(); page = 0; refresh();
          }));
          dirBtns.forEach(btn => btn.addEventListener('click', () => {
            if (evFilter === 'binding') return;
            dirBtns.forEach(b => b.classList.remove('active'));
            btn.classList.add('active'); dirFilter = btn.dataset.dir; page = 0; refresh();
          }));
        } else if (simpleBtns.length) {
          // Simple mode (tf.html, gene.html)
          simpleBtns.forEach(btn => btn.addEventListener('click', () => {
            simpleBtns.forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            const f = btn.dataset.filter;
            // Map simple filter names to ev/dir filters
            if (f === 'all') { evFilter = 'all'; dirFilter = 'all'; }
            else if (f === 'up') { evFilter = 'all'; dirFilter = 'up'; }
            else if (f === 'down') { evFilter = 'all'; dirFilter = 'down'; }
            else if (f === 'perturbation') { evFilter = 'perturbation'; dirFilter = 'all'; }
            else if (f === 'binding') { evFilter = 'binding'; dirFilter = 'all'; }
            else if (f === 'both') { evFilter = 'both'; dirFilter = 'all'; }
            else { evFilter = f; dirFilter = 'all'; }
            page = 0; refresh();
          }));
        }
      }
    }

    // Search
    if (opts.searchId) {
      const el = document.getElementById(opts.searchId);
      if (el) el.addEventListener('input', () => { searchQ = el.value.trim().toLowerCase(); page = 0; refresh(); });
    }

    refresh();
    return {
      refresh,
      getFiltered: () => applyFilter(allData),
      getFilterState: () => ({ evFilter, dirFilter, searchQ }),
    };
  }

  /**
   * Simple sortable table (no filtering).
   */
  function makeSortable(table, data, renderRow) {
    const headers = table.querySelectorAll('th.sortable');
    let sortKey = null, sortDir = 'asc';
    function render(d) {
      const tb = table.querySelector('tbody');
      tb.innerHTML = '';
      d.forEach(r => tb.appendChild(renderRow(r)));
    }
    headers.forEach(th => th.addEventListener('click', () => {
      const key = th.dataset.key, type = th.dataset.type;
      if (sortKey === key) sortDir = sortDir === 'asc' ? 'desc' : 'asc';
      else { sortKey = key; sortDir = type === 'number' ? 'desc' : 'asc'; }
      headers.forEach(h => h.classList.remove('sort-asc', 'sort-desc'));
      th.classList.add(sortDir === 'asc' ? 'sort-asc' : 'sort-desc');
      render([...data].sort((a, b) => {
        let va = a[key], vb = b[key];
        if (va == null) return 1;
        if (vb == null) return -1;
        if (type === 'number') return sortDir === 'asc' ? va - vb : vb - va;
        return sortDir === 'asc' ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
      }));
    }));
    render(data);
  }

  // ── Enrichment table ──────────────────────────────────────────────────────

  /**
   * Paginated enrichment table with source filter and search.
   */
  function makeEnrichTable(enrichTerms, opts) {
    opts = opts || {};
    const PAGE_SIZE = opts.pageSize || 10;
    const table = document.getElementById(opts.tableId || 'enrich-table');
    const emptyEl = document.getElementById(opts.emptyId || 'enrich-empty');
    const pagEl = document.getElementById(opts.pagId || 'enrich-pagination');
    const searchEl = document.getElementById(opts.searchId || 'enrich-search');
    const filterWrap = document.getElementById(opts.filterId || 'enrich-source-filters');
    if (!table) return;
    if (!enrichTerms.length) { if (emptyEl) emptyEl.style.display = ''; return; }
    table.style.display = '';

    const sources = [...new Set(enrichTerms.map(r => r.source))].sort();
    let activeSource = 'all';
    if (filterWrap) {
      filterWrap.innerHTML = '<button class="enrich-src-btn active" data-src="all">All</button>' +
        sources.map(s => '<button class="enrich-src-btn" data-src="' + s + '">' + s + '</button>').join('');
      filterWrap.addEventListener('click', e => {
        const btn = e.target.closest('.enrich-src-btn');
        if (!btn) return;
        activeSource = btn.dataset.src;
        filterWrap.querySelectorAll('.enrich-src-btn').forEach(b => b.classList.toggle('active', b === btn));
        page = 0; refresh();
      });
    }

    let sortKey = 'p_value', sortDir = 'asc', q = '', page = 0;
    function renderRow(r) {
      const tr = document.createElement('tr');
      const termLink = (r.source === 'GO:BP' || r.source === 'GO:CC' || r.source === 'GO:MF') && r.term_id
        ? '<a class="node-link" href="' + PERTURBSEQ_PREFIX + '/go/' + encodeURIComponent(r.term_id) + '">' + r.term_name + '</a>'
        : r.term_name;
      tr.innerHTML = '<td>' + termLink + '</td><td><span class="enrich-source-badge enrich-source-' + r.source.replace(':', '') + '">' + r.source + '</span></td><td>' + r.p_value.toExponential(2) + '</td><td>' + r.intersection_size + '/' + r.term_size + '</td>';
      return tr;
    }
    function refresh() {
      let d = enrichTerms;
      if (activeSource !== 'all') d = d.filter(r => r.source === activeSource);
      if (q) d = d.filter(r => r.term_name.toLowerCase().includes(q));
      d = [...d].sort((a, b) => {
        const va = a[sortKey], vb = b[sortKey];
        const type = table.querySelector('th[data-key="' + sortKey + '"]')?.dataset.type || 'string';
        if (type === 'number') return sortDir === 'asc' ? va - vb : vb - va;
        return sortDir === 'asc' ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
      });
      const totalPages = Math.max(1, Math.ceil(d.length / PAGE_SIZE));
      page = Math.min(page, totalPages - 1);
      const tb = table.querySelector('tbody');
      tb.innerHTML = '';
      d.slice(page * PAGE_SIZE, (page + 1) * PAGE_SIZE).forEach(r => tb.appendChild(renderRow(r)));
      if (pagEl) {
        pagEl.style.display = totalPages <= 1 ? 'none' : 'flex';
        pagEl.innerHTML = '<button id="enrich-prev" ' + (page === 0 ? 'disabled' : '') + '>← Prev</button><span>' + (page + 1) + ' / ' + totalPages + '</span><button id="enrich-next" ' + (page >= totalPages - 1 ? 'disabled' : '') + '>Next →</button>';
        document.getElementById('enrich-prev')?.addEventListener('click', () => { page--; refresh(); });
        document.getElementById('enrich-next')?.addEventListener('click', () => { page++; refresh(); });
      }
    }
    table.querySelectorAll('th.sortable').forEach(th => {
      th.addEventListener('click', () => {
        const key = th.dataset.key;
        if (sortKey === key) sortDir = sortDir === 'asc' ? 'desc' : 'asc';
        else { sortKey = key; sortDir = th.dataset.type === 'number' ? 'asc' : 'asc'; }
        table.querySelectorAll('th').forEach(t => t.classList.remove('sort-asc', 'sort-desc'));
        th.classList.add(sortDir === 'asc' ? 'sort-asc' : 'sort-desc');
        page = 0; refresh();
      });
    });
    if (searchEl) searchEl.addEventListener('input', () => { q = searchEl.value.trim().toLowerCase(); page = 0; refresh(); });
    table.querySelector('th[data-key="p_value"]')?.classList.add('sort-asc');
    refresh();
  }

  // ── Cytoscape filter helper ───────────────────────────────────────────────

  /**
   * Apply filter to a cytoscape instance based on edge evidence/direction.
   * Hides orphan peripheral nodes and fits the view.
   *
   * @param {Object} cy            - Cytoscape instance
   * @param {string} filter        - 'all' | 'perturbation' | 'binding' | 'both' | 'up' | 'down'
   * @param {string} [nodeType]    - Peripheral node type to check for orphans (default: auto)
   * @param {number} [padding]     - Fit padding (default: 30)
   */
  function applyCyFilter(cy, filter, nodeType, padding) {
    padding = padding != null ? padding : 30;
    cy.elements().show();
    if (filter !== 'all') {
      cy.edges().forEach(e => {
        const ev = e.data('evidence'), dir = e.data('direction');
        let hide = false;
        if (filter === 'perturbation')      hide = ev === 'binding';
        else if (filter === 'binding')      hide = ev !== 'binding' && ev !== 'both';
        else if (filter === 'both')         hide = ev !== 'both';
        else if (filter === 'up')           hide = dir !== 'up';
        else if (filter === 'down')         hide = dir !== 'down';
        if (hide) e.hide();
      });
    }
    // Hide orphan peripheral nodes
    const sel = nodeType ? 'node[type="' + nodeType + '"]' : 'node';
    cy.nodes(sel).forEach(n => {
      if (n.data('type') && !n.connectedEdges(':visible').length) n.hide();
    });
    cy.fit(cy.elements(':visible'), padding);
  }

  // ── TF status badge ───────────────────────────────────────────────────────

  /**
   * Return an HTML string for a TF/gene type badge.
   * @param {string} status  - 'Active TF' | 'TF' | 'Gene' | ''
   * @param {string} variant - 'type' (default, tf-type-* classes) | 'inline' (intro-tag classes)
   */
  function tfStatusBadge(status, variant) {
    if (variant === 'inline') {
      if (status === 'Active TF') return '<span class="intro-tag gene-tag" style="font-size:9px;padding:1px 6px;margin-left:6px">Active TF</span>';
      if (status === 'TF')        return '<span class="intro-tag term-tag" style="font-size:9px;padding:1px 6px;margin-left:6px">TF</span>';
      return '';
    }
    if (status === 'Active TF') return '<span class="tf-type-active">Active TF</span>';
    if (status === 'TF')        return '<span class="tf-type-tf">TF</span>';
    return '<span style="color:var(--text-muted);font-size:11px">Gene</span>';
  }

  // ── Public API ────────────────────────────────────────────────────────────
  return {
    MODULE_COLORS,
    EDGE_STYLES_FULL,
    EDGE_STYLES_GC,
    hubNodeStyle,
    peripheralNodeStyle,
    placeSector,
    initNetwork,
    makeFilterSortTable,
    makeSortable,
    makeEnrichTable,
    applyCyFilter,
    tfStatusBadge,
  };
})();
