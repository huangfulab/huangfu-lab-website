/**
 * drawExprChart — lightweight TPM expression line chart.
 * Shared by tf_gene_link.html and any other template needing side-by-side mini charts.
 * Requires D3 v7 to be loaded first.
 *
 * @param {string} svgId  - ID of the <svg> element to draw into
 * @param {Array}  data   - Array of {timepoint, mean_tpm}
 * @param {string} color  - Stroke/fill colour
 */
function drawExprChart(svgId, data, color) {
  const svgEl = document.getElementById(svgId);
  if (!svgEl || !data.length) return;

  const margin = { top: 12, right: 10, bottom: 38, left: 44 };
  const W  = svgEl.parentElement.clientWidth - 2;
  const H  = 160;
  const iW = W - margin.left - margin.right;
  const iH = H - margin.top  - margin.bottom;

  svgEl.setAttribute('width', W);
  svgEl.setAttribute('height', H);

  const tps  = data.map(d => d.timepoint);
  const vals = data.map(d => d.mean_tpm || 0);
  const xMax = Math.max(...vals);
  const pad  = xMax * 0.08 || 0.2;

  const x = d3.scalePoint().domain(tps).range([0, iW]).padding(0.3);
  const y = d3.scaleLinear().domain([0, xMax + pad]).range([iH, 0]).nice();

  const svg = d3.select(svgEl).append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`);

  svg.append('g').attr('class', 'grid')
    .call(d3.axisLeft(y).tickSize(-iW).tickFormat(''));
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${iH})`)
    .call(d3.axisBottom(x).tickFormat(t => t.replace('_', ' ')))
    .selectAll('text').attr('transform', 'rotate(-35)').style('text-anchor', 'end');
  svg.append('g').attr('class', 'axis')
    .call(d3.axisLeft(y).ticks(4));
  svg.append('text')
    .attr('transform', 'rotate(-90)').attr('x', -iH / 2).attr('y', -32)
    .attr('text-anchor', 'middle').attr('font-size', 9).attr('fill', 'var(--text-muted)')
    .text('TPM');

  svg.append('path').datum(data).attr('class', 'area-fill')
    .attr('fill', color)
    .attr('d', d3.area()
      .x(d => x(d.timepoint))
      .y0(iH).y1(d => y(d.mean_tpm || 0))
      .curve(d3.curveCatmullRom.alpha(0.5)));

  svg.append('path').datum(data).attr('class', 'mean-line')
    .attr('stroke', color)
    .attr('d', d3.line()
      .x(d => x(d.timepoint))
      .y(d => y(d.mean_tpm || 0))
      .curve(d3.curveCatmullRom.alpha(0.5)));

  // Shared tooltip — one div per page, reused across all drawExprChart calls
  let tip = d3.select('#_expr-chart-tip');
  if (tip.empty()) {
    tip = d3.select('body').append('div').attr('id', '_expr-chart-tip')
      .style('position', 'fixed').style('pointer-events', 'none')
      .style('background', 'var(--bg-elevated)').style('border', '1px solid var(--border)')
      .style('border-radius', '6px').style('padding', '5px 9px')
      .style('font-size', '11px').style('color', 'var(--text)')
      .style('box-shadow', '0 2px 8px rgba(0,0,0,0.15)')
      .style('display', 'none').style('z-index', '9999').style('white-space', 'nowrap');
  }

  data.forEach(d => {
    svg.append('circle')
      .attr('cx', x(d.timepoint)).attr('cy', y(d.mean_tpm || 0)).attr('r', 4.5)
      .attr('fill', color).attr('stroke', '#fff').attr('stroke-width', 1.5)
      .style('cursor', 'default')
      .on('mouseover', () => {
        tip.html(`<strong>${d.timepoint.replace('_', ' ')}</strong><br>${(d.mean_tpm || 0).toFixed(2)} TPM`)
           .style('display', 'block');
      })
      .on('mousemove', ev => tip.style('left', (ev.clientX + 14) + 'px').style('top', (ev.clientY - 44) + 'px'))
      .on('mouseout',  () => tip.style('display', 'none'));
  });
}
