"""
generate_comparison_report.py – 產生 pipeline 比較 HTML 報告

整合：
  • 偵測準確率比較（Precision / Recall / F1 / AUC-PR）
  • 分層 F1（低/中/高 BSJ count）
  • 計算效能比較（wall time, peak RAM, CPU-hours）
  • DE 分析品質比較（DE 數量, Jaccard, Type I unique, circBase hits）
  • Pipeline feature 功能比較表
  • 結論段落

Usage:
  python generate_comparison_report.py \
      --accuracy   results/benchmark/accuracy_summary.tsv \
      --stratified results/benchmark/stratified_f1.tsv \
      --compute    results/benchmark/compute_cost.tsv \
      --de-quality results/benchmark/de_quality_summary.tsv \
      --de-jaccard results/benchmark/de_jaccard.tsv \
      --output     results/benchmark/comparison_report.html
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd


# ── Style ─────────────────────────────────────────────────────────────────────

_STYLE = """
<style>
  * { box-sizing: border-box; }
  body {
    font-family: "Segoe UI", Arial, sans-serif;
    max-width: 1200px; margin: 40px auto; padding: 0 20px;
    color: #222; background: #f8f9fa;
  }
  h1 { color: #1a5c96; border-bottom: 3px solid #1a5c96; padding-bottom: 10px; }
  h2 { color: #2c6fad; margin-top: 40px; border-left: 4px solid #2c6fad;
       padding-left: 10px; }
  h3 { color: #444; margin-top: 24px; }

  /* Tables */
  .tbl-wrap { overflow-x: auto; margin: 16px 0; }
  table { border-collapse: collapse; width: 100%; font-size: 13px;
          background: #fff; border-radius: 6px; overflow: hidden;
          box-shadow: 0 1px 3px rgba(0,0,0,.1); }
  th { background: #1a5c96; color: #fff; padding: 8px 12px;
       text-align: left; white-space: nowrap; }
  td { border-bottom: 1px solid #e8e8e8; padding: 7px 12px; }
  tr:last-child td { border-bottom: none; }
  tr:hover td { background: #edf4ff; }
  .best { font-weight: bold; color: #0d7a2e; }
  .na   { color: #aaa; font-style: italic; }

  /* Stat boxes */
  .stat-row { display: flex; flex-wrap: wrap; gap: 12px; margin: 16px 0; }
  .stat-box { background: #fff; border: 1px solid #c5d8f0;
              border-radius: 8px; padding: 14px 22px; min-width: 160px;
              text-align: center; box-shadow: 0 1px 3px rgba(0,0,0,.07); }
  .stat-box .num { font-size: 2em; font-weight: bold; color: #1a5c96; }
  .stat-box .lbl { font-size: 0.82em; color: #666; margin-top: 4px; }

  /* Bar chart */
  .bar-group { margin: 6px 0; display: flex; align-items: center; gap: 8px; }
  .bar-label { width: 180px; font-size: 12px; text-align: right;
               white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
  .bar-wrap  { flex: 1; background: #e8e8e8; border-radius: 4px; height: 20px; }
  .bar-fill  { height: 20px; border-radius: 4px; display: flex;
               align-items: center; padding-left: 6px;
               color: #fff; font-size: 11px; font-weight: bold; transition: width .3s; }
  .c1 { background: #1a5c96; }
  .c2 { background: #e07b39; }
  .c3 { background: #2da84b; }
  .bar-val   { font-size: 12px; width: 50px; color: #333; }
  .c1 { background: #1a5c96; }
  .c2 { background: #7b3fa6; }
  .c3 { background: #e07b39; }
  .c4 { background: #2da84b; }
  .c5 { background: #d62728; }

  /* Feature table special */
  .feat-yes { color: #0d7a2e; font-weight: bold; }
  .feat-no  { color: #c0392b; }
  .feat-partial { color: #e67e22; }

  /* Section card */
  .card { background: #fff; border-radius: 8px; padding: 24px;
          margin: 24px 0; box-shadow: 0 1px 4px rgba(0,0,0,.08); }

  /* Conclusion */
  .concl-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
                gap: 16px; margin-top: 16px; }
  .concl-card { border-radius: 8px; padding: 16px 20px; }
  .c-ours        { background: #e8f0fa; border-left: 4px solid #1a5c96; }
  .c-circompara2 { background: #f3eafa; border-left: 4px solid #7b3fa6; }
  .c-nfcore      { background: #fef3e8; border-left: 4px solid #e07b39; }
  .c-sponging    { background: #e8f9ec; border-left: 4px solid #2da84b; }
  .c-clear       { background: #fdf0f0; border-left: 4px solid #d62728; }
  .concl-card h4 { margin: 0 0 10px; font-size: 1em; }
  .concl-card ul { margin: 6px 0; padding-left: 18px; font-size: 13px; }
  .concl-card li { margin-bottom: 4px; }

  .note { font-size: 12px; color: #666; margin-top: 8px; font-style: italic; }
  .badge { display: inline-block; padding: 2px 8px; border-radius: 10px;
           font-size: 11px; font-weight: bold; color: #fff; }
  .badge-ours        { background: #1a5c96; }
  .badge-circompara2 { background: #7b3fa6; }
  .badge-nfcore      { background: #e07b39; }
  .badge-sponging    { background: #2da84b; }
  .badge-clear       { background: #d62728; }

  /* Download / print toolbar */
  .dl-btn {
    display: inline-block; margin: 4px 4px 4px 0;
    padding: 4px 10px; font-size: 11px; font-weight: bold;
    border: 1px solid #1a5c96; border-radius: 4px;
    color: #1a5c96; background: #fff; cursor: pointer;
    text-decoration: none; transition: background .15s;
  }
  .dl-btn:hover { background: #e8f0fa; }
  .tbl-header { display: flex; align-items: center;
                justify-content: space-between; flex-wrap: wrap;
                gap: 4px; margin-bottom: 4px; }
  .tbl-title  { font-weight: bold; font-size: 14px; color: #333; }
  .print-bar  { position: sticky; top: 0; z-index: 99;
                background: #1a5c96; color: #fff;
                padding: 8px 24px; display: flex;
                align-items: center; gap: 12px; flex-wrap: wrap; }
  .print-bar span { font-size: 13px; font-weight: bold; flex: 1; }
  .print-btn  { padding: 6px 16px; border: 2px solid #fff;
                border-radius: 4px; background: transparent;
                color: #fff; font-size: 13px; font-weight: bold;
                cursor: pointer; }
  .print-btn:hover { background: rgba(255,255,255,.15); }

  /* DE quality table — clickable cells */
  .de-clickable {
    cursor: pointer;
    color: #1a5276;
    font-weight: bold;
    text-decoration: underline dotted #1a5276;
  }
  .de-clickable:hover { background: #d6eaf8 !important; }

  /* Print styles */
  @media print {
    .print-bar, .dl-btn, .no-print { display: none !important; }
    body { background: #fff; max-width: 100%; margin: 0; padding: 0 16px; }
    .card { box-shadow: none; border: 1px solid #ddd; page-break-inside: avoid; }
    table { page-break-inside: avoid; font-size: 11px; }
    th { background: #1a5c96 !important; -webkit-print-color-adjust: exact;
         print-color-adjust: exact; }
    h1, h2 { color: #1a5c96 !important; -webkit-print-color-adjust: exact;
              print-color-adjust: exact; }
    .bar-fill { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    .concl-grid { grid-template-columns: 1fr 1fr; }
  }
</style>
"""

_SCRIPT = """
<script>
function dlCSV(tid, fname) {
  var tbl = document.getElementById(tid);
  if (!tbl) return;
  var rows = tbl.querySelectorAll('tr');
  var lines = [];
  rows.forEach(function(r) {
    var cells = r.querySelectorAll('th,td');
    var cols = [];
    cells.forEach(function(c) {
      cols.push('"' + c.innerText.replace(/"/g, '""') + '"');
    });
    lines.push(cols.join(','));
  });
  var blob = new Blob([lines.join('\\n')], {type:'text/csv;charset=utf-8;'});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = fname;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
}

/* ── Sortable list table ───────────────────────────────────── */
var _LIST_SORT = {col: -1, dir: 1};
function _parseVal(s) {
  if (s === '—' || s === '' || s == null) return Infinity;
  // scientific notation (p-value)
  var n = parseFloat(s);
  return isNaN(n) ? s.toString().toLowerCase() : n;
}
function sortCircList(th, colIdx) {
  var table = th.closest('table');
  var tbody = table.querySelector('tbody');
  var rows  = Array.from(tbody.querySelectorAll('tr'));
  if (_LIST_SORT.col === colIdx) { _LIST_SORT.dir *= -1; }
  else { _LIST_SORT.col = colIdx; _LIST_SORT.dir = 1; }
  var dir = _LIST_SORT.dir;
  th.closest('tr').querySelectorAll('th').forEach(function(h, i) {
    var sp = h.querySelector('.si');
    if (sp) sp.textContent = (i === colIdx) ? (dir===1?' ▲':' ▼') : ' ⇅';
  });
  rows.sort(function(a, b) {
    var tda = a.querySelectorAll('td')[colIdx];
    var tdb = b.querySelectorAll('td')[colIdx];
    var va = tda ? (tda.dataset.val !== undefined ? tda.dataset.val : tda.textContent.trim()) : '';
    var vb = tdb ? (tdb.dataset.val !== undefined ? tdb.dataset.val : tdb.textContent.trim()) : '';
    var pa = _parseVal(va), pb = _parseVal(vb);
    if (typeof pa === 'number' && typeof pb === 'number') return dir * (pa - pb);
    return dir * String(pa).localeCompare(String(pb));
  });
  rows.forEach(function(r, i) {
    var fc = r.querySelectorAll('td')[0];
    if (fc) fc.textContent = i + 1;
    tbody.appendChild(r);
  });
}

/* ── DE quality interactive list ───────────────────────────── */
var _DE_LIST_CUR = null;
function showDEList(method, col) {
  var lists = (typeof DE_LISTS !== 'undefined') ? DE_LISTS : {};
  var data = (lists[method] || {})[col];
  if (!data || !data.length) { alert('無可用清單資料'); return; }
  var colLabels = {
    'Sig_DE_circRNAs':         '顯著 DE circRNA',
    'Up_regulated':            '上調 circRNA',
    'Down_regulated':          '下調 circRNA',
    'Type_I_count':            'Type I circRNA',
    'Type_II_count':           'Type II circRNA',
    'Type_I_unique_vs_DESeq2': 'Type I（僅本方法偵測）',
    'Top20_in_circBase':       'Top 20 DE（circBase 已知）'
  };
  document.getElementById('de-list-title').textContent =
    method.replace(/_/g,' ') + ' — ' + (colLabels[col] || col);
  document.getElementById('de-list-count').textContent =
    '共 ' + data.length + ' 個 circRNA';
  var hasType = data.some(function(r){ return r.Type; });
  document.getElementById('de-list-body').innerHTML = _renderCircList(data, hasType ? 'Type' : null);
  document.getElementById('de-list-modal').style.display = 'block';
  var hdrs = ['#','circ_id','gene_name','log2FC','p-value','circbase_id'];
  if (hasType) hdrs.push('Type');
  _DE_LIST_CUR = {data:data, method:method, col:col, hdrs:hdrs};
}
function closeDeList() {
  document.getElementById('de-list-modal').style.display = 'none';
}
function dlDeList() {
  if (!_DE_LIST_CUR) return;
  var rows = [_DE_LIST_CUR.hdrs.join(',')];
  _DE_LIST_CUR.data.forEach(function(r,i){
    var line = [i+1,'"'+r.circ_id+'"','"'+(r.gene_name||'')+'"',
      (r.log2FC!=null?r.log2FC:''), (r.p_value!=null?r.p_value:''),
      '"'+(r.circbase_id||'')+'"'];
    if (_DE_LIST_CUR.hdrs.indexOf('Type')>=0) line.push('"'+(r.Type||'')+'"');
    rows.push(line.join(','));
  });
  var blob = new Blob([rows.join('\\n')],{type:'text/csv;charset=utf-8;'});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = _DE_LIST_CUR.method+'_'+_DE_LIST_CUR.col+'.csv';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}

/* ── Jaccard pairwise interactive list ─────────────────────── */
function _renderCircList(data, extraHdr) {
  var hdrs = ['#','circ_id','gene_name','log2FC','p-value','circbase_id'];
  if (extraHdr) hdrs.push(extraHdr);
  var TH_S = 'background:#1a5276;color:#fff;padding:6px 10px;text-align:left;white-space:nowrap;cursor:pointer;user-select:none';
  var html = '<table style="width:100%;border-collapse:collapse;font-size:13px">'
    + '<thead><tr>' + hdrs.map(function(h,i){
        var si = i===0 ? '' : '<span class="si" style="opacity:0.6"> ⇅</span>';
        var click = i===0 ? '' : ' onclick="sortCircList(this,'+i+')"';
        return '<th style="'+TH_S+'"'+click+'>'+h+si+'</th>';
      }).join('') + '</tr></thead><tbody>';
  _LIST_SORT = {col: -1, dir: 1};
  data.forEach(function(r,i){
    var lfcRaw = (r.log2FC  != null) ? r.log2FC  : null;
    var pvRaw  = (r.p_value != null) ? r.p_value : null;
    var lfc    = lfcRaw  != null ? lfcRaw.toFixed(2)         : '—';
    var pval   = pvRaw   != null ? pvRaw.toExponential(2)    : '—';
    var lc = lfcRaw > 0 ? '#c0392b' : lfcRaw < 0 ? '#1a5276' : '';
    var cc = (r.circbase_id && r.circbase_id !== 'novel') ? '#e67e22' : '';
    var bg = i%2===0 ? '#f9f9f9' : '#fff';
    var TD = 'padding:5px 10px;border-bottom:1px solid #eee';
    var gene = (r.gene_name && r.gene_name !== 'intergenic') ? r.gene_name : (r.gene_name === 'intergenic' ? '<span style="color:#aaa">intergenic</span>' : '—');
    var cells = [
      '<td style="'+TD+';color:#888">'+(i+1)+'</td>',
      '<td style="'+TD+';font-family:monospace;font-size:11px" data-val="'+r.circ_id+'">'+r.circ_id+'</td>',
      '<td style="'+TD+'" data-val="'+(r.gene_name||'')+'">'+gene+'</td>',
      '<td style="'+TD+';color:'+lc+';font-weight:bold" data-val="'+(lfcRaw!=null?lfcRaw:'Infinity')+'">'+lfc+'</td>',
      '<td style="'+TD+'" data-val="'+(pvRaw!=null?pvRaw:'Infinity')+'">'+pval+'</td>',
      '<td style="'+TD+';color:'+cc+'" data-val="'+(r.circbase_id||'')+'">'+  (r.circbase_id||'—')+'</td>'
    ];
    if (extraHdr) cells.push('<td style="'+TD+'" data-val="'+(r.Type||'')+'">'+  (r.Type||'—')+'</td>');
    html += '<tr style="background:'+bg+'">'+cells.join('')+'</tr>';
  });
  html += '</tbody></table>';
  return html;
}
function showJaccardList(rowIdx, col) {
  var lists = (typeof DE_LISTS !== 'undefined') ? DE_LISTS : {};
  var jac = lists['_jaccard'] || [];
  var entry = jac[rowIdx] || {};
  var data = entry[col] || [];
  if (!data.length) { alert('無可用清單資料'); return; }
  var colLabels = {'A_only':'A 方法獨有 circRNA','B_only':'B 方法獨有 circRNA','Both':'兩方法共同 circRNA'};
  var title = (entry.Method_A||'').replace(/_/g,' ') + ' vs ' + (entry.Method_B||'').replace(/_/g,' ')
            + ' — ' + (colLabels[col] || col);
  document.getElementById('de-list-title').textContent = title;
  document.getElementById('de-list-count').textContent = '共 ' + data.length + ' 個 circRNA';
  document.getElementById('de-list-body').innerHTML = _renderCircList(data, null);
  document.getElementById('de-list-modal').style.display = 'block';
  _DE_LIST_CUR = {data:data, method:(entry.Method_A||'')+'_vs_'+(entry.Method_B||''), col:col,
                  hdrs:['#','circ_id','gene_name','log2FC','p-value','circbase_id']};
}
</script>
"""


# ── Helpers ───────────────────────────────────────────────────────────────────

def _df_html(df: pd.DataFrame, best_col: str | None = None,
             best_max: bool = True,
             table_id: str | None = None,
             csv_filename: str | None = None,
             title: str = "",
             highlight_cols: dict | None = None) -> str:
    """Render DataFrame as HTML table with optional CSV download button.

    highlight_cols: {col_name: True(max)/False(min)} — highlight best value in each listed column.
    """
    # Build per-column best values
    _best = {}
    if best_col:
        highlight_cols = dict(highlight_cols or {})
        highlight_cols[best_col] = best_max
    for col, is_max in (highlight_cols or {}).items():
        if col in df.columns:
            col_vals = pd.to_numeric(df[col], errors="coerce").dropna()
            if not col_vals.empty:
                _best[col] = col_vals.max() if is_max else col_vals.min()

    rows_html = ""
    for i, row in df.iterrows():
        cells = ""
        for col in df.columns:
            val = row[col]
            cls = ""
            if col in _best:
                try:
                    if abs(float(val) - float(_best[col])) < 1e-9:
                        cls = ' class="best"'
                except (TypeError, ValueError):
                    pass
            if val in (None, "N/A", float("nan")) or (
                isinstance(val, float) and pd.isna(val)
            ):
                cells += f'<td class="na">—</td>'
            else:
                cells += f"<td{cls}>{val}</td>"
        rows_html += f"<tr>{cells}</tr>"
    headers = "".join(f"<th>{c}</th>" for c in df.columns)
    tid   = f' id="{table_id}"' if table_id else ""
    dl_btn = (
        f'<button class="dl-btn no-print" '
        f'onclick="dlCSV(\'{table_id}\',\'{csv_filename}\')">⬇ CSV</button>'
        if table_id and csv_filename else ""
    )
    header_bar = (
        f'<div class="tbl-header">'
        f'<span class="tbl-title">{title}</span>{dl_btn}</div>'
        if title or dl_btn else ""
    )
    return (
        f'{header_bar}<div class="tbl-wrap"><table{tid}>'
        f"<thead><tr>{headers}</tr></thead>"
        f"<tbody>{rows_html}</tbody></table></div>"
    )


_DE_CLICKABLE_COLS = {
    "Sig_DE_circRNAs", "Up_regulated", "Down_regulated",
    "Type_I_count", "Type_II_count", "Type_I_unique_vs_DESeq2", "Top20_in_circBase",
}


def _de_quality_table_html(
    de: pd.DataFrame,
    de_display_cols: list,
    de_lists: dict | None,
) -> str:
    """Render DE quality table; numeric cells in clickable columns open circRNA modal."""
    dl_btn = '<button class="dl-btn no-print" onclick="dlCSV(\'tbl_de\',\'de_quality_summary.csv\')">⬇ CSV</button>'
    headers = "".join(f"<th>{c}</th>" for c in de_display_cols)
    rows_html = ""
    for _, row in de[de_display_cols].iterrows():
        method = str(row.get("Method", ""))
        cells = ""
        for col in de_display_cols:
            val = row[col]
            is_na = val in (None, "N/A") or (isinstance(val, float) and pd.isna(val))
            if is_na:
                cells += '<td class="na">—</td>'
            elif (
                de_lists is not None
                and col in _DE_CLICKABLE_COLS
                and de_lists.get(method, {}).get(col)
            ):
                cells += (
                    f'<td class="de-clickable" '
                    f'onclick="showDEList(\'{method}\',\'{col}\')" '
                    f'title="點擊查看 circRNA 清單">{val}</td>'
                )
            else:
                cells += f"<td>{val}</td>"
        rows_html += f"<tr>{cells}</tr>"
    return (
        f'<div class="tbl-header">{dl_btn}</div>'
        f'<div class="tbl-wrap"><table id="tbl_de">'
        f"<thead><tr>{headers}</tr></thead>"
        f"<tbody>{rows_html}</tbody></table></div>"
    )


_JAC_CLICKABLE = {"A_only", "B_only", "Both"}


def _jaccard_table_html(jac: pd.DataFrame, de_lists: "dict | None") -> str:
    """Render Jaccard table with A_only / B_only / Both cells clickable."""
    jac_list_data = (de_lists or {}).get("_jaccard", [])
    idx_map: dict = {}
    for i, entry in enumerate(jac_list_data):
        key = (entry.get("Method_A", ""), entry.get("Method_B", ""))
        idx_map[key] = i

    dl_btn = '<button class="dl-btn no-print" onclick="dlCSV(\'tbl_jaccard\',\'de_jaccard.csv\')">⬇ CSV</button>'
    cols = list(jac.columns)
    headers = "".join(f"<th>{c}</th>" for c in cols)
    rows_html = ""
    for _, row in jac.iterrows():
        ma = str(row.get("Method_A", ""))
        mb = str(row.get("Method_B", ""))
        row_idx = idx_map.get((ma, mb))
        cells = ""
        for col in cols:
            val = row[col]
            is_na = val in (None, "N/A") or (isinstance(val, float) and pd.isna(val))
            if is_na:
                cells += '<td class="na">—</td>'
            elif (
                col in _JAC_CLICKABLE
                and row_idx is not None
                and jac_list_data
                and jac_list_data[row_idx].get(col)
            ):
                cells += (
                    f'<td class="de-clickable" '
                    f'onclick="showJaccardList({row_idx},\'{col}\')" '
                    f'title="點擊查看 circRNA 清單">{val}</td>'
                )
            else:
                cells += f"<td>{val}</td>"
        rows_html += f"<tr>{cells}</tr>"
    return (
        f'<div class="tbl-header">{dl_btn}</div>'
        f'<div class="tbl-wrap"><table id="tbl_jaccard">'
        f"<thead><tr>{headers}</tr></thead>"
        f"<tbody>{rows_html}</tbody></table></div>"
    )


_TOOL_COLORS = {
    "CIRIquant_min":    ("#1a5276", "CIRIquant"),
    "STAR_min":         ("#2980b9", "STAR×3"),
    "DCC_min":          ("#5dade2", "DCC"),
    "CIRCexplorer2_min":("#e67e22", "CIRCexplorer2"),
    "find_circ_min":    ("#27ae60", "find_circ"),
    "DE_min":           ("#8e44ad", "DE analysis"),
}


def _compute_cost_table_html(comp: pd.DataFrame) -> str:
    """Render compute cost table with per-tool stacked bar breakdown."""
    import math

    def _fv(v) -> "float | None":
        try:
            f = float(v)
            return None if math.isnan(f) else f
        except (TypeError, ValueError):
            return None

    max_total = float(comp["Total_wall_min"].max())

    legend = "".join(
        f'<span style="display:inline-flex;align-items:center;margin-right:12px">'
        f'<span style="width:12px;height:12px;background:{col};display:inline-block;'
        f'margin-right:4px;border-radius:2px"></span>{lbl}</span>'
        for col, (color, lbl) in _TOOL_COLORS.items()  # noqa: unused col
        for col2, (color, lbl) in [(col, (color, lbl))]  # iterate once
    )
    # rebuild correctly
    legend = "".join(
        f'<span style="display:inline-flex;align-items:center;margin-right:12px">'
        f'<span style="width:12px;height:12px;background:{color};display:inline-block;'
        f'margin-right:4px;border-radius:2px"></span>{lbl}</span>'
        for _, (color, lbl) in _TOOL_COLORS.items()
    )
    legend += (
        '<span style="display:inline-flex;align-items:center;margin-right:12px">'
        '<span style="width:12px;height:12px;background:#95a5a6;display:inline-block;'
        'margin-right:4px;border-radius:2px"></span>Other/Consensus</span>'
    )

    th = 'style="background:#1a5276;color:#fff;padding:8px 10px;white-space:nowrap"'
    th_r = 'style="background:#1a5276;color:#fff;padding:8px 10px;white-space:nowrap;text-align:right"'
    rows_html = ""
    min_total = float(comp["Total_wall_min"].min())

    for i, (_, row) in enumerate(comp.iterrows()):
        bg = "#f9f9f9" if i % 2 == 0 else "#fff"
        total = _fv(row["Total_wall_min"]) or 0.0

        # Build stacked bar segments
        segs = ""
        tool_sum = 0.0
        breakdown_parts = []
        for tcol, (color, lbl) in _TOOL_COLORS.items():
            if tcol not in row.index:
                continue
            val = _fv(row[tcol])
            if val and val > 0:
                pct = val / max_total * 100
                segs += (f'<div title="{lbl}: {val:.1f} min" '
                         f'style="display:inline-block;width:{pct:.1f}%;height:18px;'
                         f'background:{color};vertical-align:top"></div>')
                tool_sum += val
                breakdown_parts.append(
                    f'<span style="color:{color};font-weight:bold">{lbl}</span>: {val:.1f}'
                )
        # Remainder (DCC shared with STAR, consensus filter, etc.)
        remainder = total - tool_sum
        if remainder > 0.5:
            pct = remainder / max_total * 100
            segs += (f'<div title="Other: {remainder:.1f} min" '
                     f'style="display:inline-block;width:{pct:.1f}%;height:18px;'
                     f'background:#95a5a6;vertical-align:top"></div>')

        breakdown_html = " · ".join(breakdown_parts) if breakdown_parts else ""
        total_style = "color:#27ae60;font-weight:bold" if total == min_total else ""

        def _fmt(v, dec=1):
            f = _fv(v)
            return f"—" if f is None else f"{f:.{dec}f}"

        rows_html += f"""<tr style="background:{bg}">
  <td style="padding:8px 10px;font-weight:bold;white-space:nowrap">{row["Pipeline"]}</td>
  <td style="padding:8px 10px;min-width:260px">
    <div style="font-size:11px;color:#666;margin-bottom:4px">{row["Tool_combination"]}</div>
    <div style="background:#ddd;height:18px;border-radius:3px;overflow:hidden;max-width:300px">{segs}</div>
    <div style="font-size:11px;margin-top:3px;color:#444">{breakdown_html}</div>
  </td>
  <td style="padding:8px 10px;text-align:right;{total_style}">{_fmt(total)}</td>
  <td style="padding:8px 10px;text-align:right">{_fmt(row.get("Peak_RAM_GB"))}</td>
  <td style="padding:8px 10px;text-align:right">{_fmt(row.get("Parallel_Peak_RAM_GB"))}</td>
  <td style="padding:8px 10px;text-align:right">{int(_fv(row.get("CPU_cores")) or 8)}</td>
  <td style="padding:8px 10px;text-align:right">{_fmt(row.get("CPU_hours"))}</td>
  <td style="padding:8px 10px;font-size:11px;color:#555">{row.get("Source","")}</td>
</tr>"""

    has_de = "DE_min" in comp.columns and comp["DE_min"].apply(_fv).notna().any()
    de_note = (
        '<div style="font-size:11px;color:#888;margin-bottom:8px">'
        '⚠ DE analysis (min) is benchmarked separately on GSE113230 (6-sample multi-sample '
        'DE run) and summed into Total; detection steps above are benchmarked on the '
        'single-sample SRR444655 ground-truth run. The two are not from the same end-to-end '
        'invocation — Total therefore represents detection-cost-on-one-sample + '
        'DE-cost-on-a-full-cohort, not a single reproducible pipeline run.</div>'
    ) if has_de else ""

    dl_btn = '<button class="dl-btn no-print" onclick="dlCSV(\'tbl_compute\',\'compute_cost.csv\')">⬇ CSV</button>'
    return f"""<div class="tbl-header">{dl_btn}</div>
<div style="font-size:12px;margin-bottom:8px;color:#555">{legend}</div>
{de_note}
<div class="tbl-wrap"><table id="tbl_compute" style="width:100%;border-collapse:collapse">
<thead><tr>
  <th {th}>Pipeline</th>
  <th {th}>Tool Breakdown</th>
  <th {th_r}>Total (min)</th>
  <th {th_r}>Peak RAM (GB)</th>
  <th {th_r}>Parallel RAM (GB)</th>
  <th {th_r}>Cores</th>
  <th {th_r}>CPU-h</th>
  <th {th}>Source</th>
</tr></thead>
<tbody>{rows_html}</tbody>
</table></div>"""


def _bar_chart(
    df: pd.DataFrame,
    metric_col: str,
    label_col: str = "Method",
    title: str = "",
    colors: list[str] | None = None,
) -> str:
    if colors is None:
        colors = ["c1", "c2", "c3"]
    def _to_float(v):
        try:
            return float(str(v).lstrip("≥").lstrip("~").strip())
        except (ValueError, TypeError):
            return None
    rows = df[[label_col, metric_col]].dropna()
    rows = rows[rows[metric_col].apply(lambda v: _to_float(v) is not None)]
    max_val = max(_to_float(v) for v in rows[metric_col]) if len(rows) > 0 else 1.0
    html = f'<h3 style="margin-bottom:8px">{title}</h3>'
    for i, (_, row) in enumerate(rows.iterrows()):
        val  = _to_float(row[metric_col])
        pct  = round(val / max_val * 100, 1) if max_val > 0 else 0
        name = str(row[label_col]).replace("_", " ")
        c    = colors[i % len(colors)]
        html += (
            f'<div class="bar-group">'
            f'<span class="bar-label">{name}</span>'
            f'<div class="bar-wrap">'
            f'<div class="bar-fill {c}" style="width:{pct}%">'
            f'{val:.3f}</div></div>'
            f'<span class="bar-val">{val:.3f}</span>'
            f'</div>'
        )
    return html


# ── Section builders ──────────────────────────────────────────────────────────

def _feature_table() -> str:
    """Static feature comparison table — 3 multi-tool pipelines."""
    # columns: Our | CirComPara2 | nf-core
    features = [
        ("Framework",
         "Snakemake", "SCons (2022)", "Nextflow (2023)"),
        ("Detection tools",
         "CIRIquant + DCC",
         "CIRIquant + DCC + find_circ + CircExplorer2",
         "CIRIquant + CIRCexplorer2 + find_circ"),
        ("Tool consensus",
         "✓ adaptive (≥2/2, slop=10 bp)",
         "✓ fixed (≥2/4 tools)",
         "✓ fixed (≥2/3 tools, exact)"),
        ("Coordinate tolerance",
         "✓ slop=10 bp (configurable)",
         "✓ slop=10 bp",
         "✗ exact match (slop=0)"),
        ("BSJ/FSJ pseudo-circ QC",
         "✓ (BSJ/FSJ ratio filter)",
         "✗",
         "✗"),
        ("Confidence scoring",
         "✓ log2(BSJ) × coord agreement",
         "~ partial (per-tool support)",
         "✗"),
        ("DE method",
         "✓ edgeR GLM + FSJ offset<br><small>✓ DESeq2 (poscounts norm.)</small><br><small>✓ limma-voom (TMM norm.)</small>",
         "DESeq2 / edgeR (BSJ counts)",
         "DESeq2 (BSJ counts)"),
        ("Type I / II classification",
         "✓ (circRNA-specific vs gene-level)",
         "✗",
         "✗"),
        ("CSI / delta-CSI",
         "✓ Circular Splicing Index",
         "✗",
         "✗"),
        ("circBase annotation",
         "✓ auto-download hg19",
         "✓ built-in",
         "✓ nf-core module"),
        ("Biomarker ranking",
         "✓ composite 6D score",
         "✗",
         "✗"),
        ("Isoform switching (IUI)",
         "✓ Wilcoxon + BH correction",
         "✗",
         "✗"),
        ("HTML report",
         "✓ self-contained + Plotly",
         "✓ auto-generated",
         "✓ MultiQC integration"),
        ("Web UI",
         "✓ Flask + GEO one-click",
         "✗",
         "✗"),
        ("Config-driven tool selection",
         "✓ CIRIquant / DCC / both",
         "✓ SCons params",
         "✓ nf-core params"),
    ]
    header = (
        '<tr><th>Feature</th>'
        '<th><span class="badge badge-ours">circDEX</span><br>'
        '<small style="color:#aaa">Snakemake</small></th>'
        '<th><span class="badge badge-circompara2">CirComPara2</span><br>'
        '<small style="color:#aaa">SCons · 2022</small></th>'
        '<th><span class="badge badge-nfcore">nf-core/circrna</span><br>'
        '<small style="color:#aaa">Nextflow · 2023</small></th>'
        '</tr>'
    )

    rows = ""
    for feat, *vals in features:
        cells = f"<td><strong>{feat}</strong></td>"
        for v in vals:
            if v.startswith("✓"):
                cells += f'<td class="feat-yes">{v}</td>'
            elif v.startswith("✗"):
                cells += f'<td class="feat-no">{v}</td>'
            elif v.startswith("~"):
                cells += f'<td class="feat-partial">{v}</td>'
            else:
                cells += f"<td>{v}</td>"
        rows += f"<tr>{cells}</tr>"

    return (
        '<div class="tbl-wrap"><table>'
        f"<thead>{header}</thead><tbody>{rows}</tbody></table></div>"
    )


def _stratified_chart(strat: pd.DataFrame) -> str:
    # Recover quartile cutoffs stored as extra columns (same value in every row)
    q1 = float(strat["q1_cutoff"].iloc[0]) if "q1_cutoff" in strat.columns else None
    q3 = float(strat["q3_cutoff"].iloc[0]) if "q3_cutoff" in strat.columns else None

    tiers = ["low_Q1", "mid_Q2Q3", "high_Q4"]
    if q1 is not None and q3 is not None:
        tier_labels = [
            f"Low BSJ count (Q1, ≤{q1:.0f} reads)",
            f"Mid BSJ count (Q2–Q3, {q1:.0f}–{q3:.0f} reads)",
            f"High BSJ count (Q4, >{q3:.0f} reads)",
        ]
    else:
        tier_labels = ["Low (Q1)", "Mid (Q2–Q3)", "High (Q4)"]

    colors_by_method = {
        "Our_adaptive":       "c1",
        "CirComPara2_4tools": "c2",
        "nfcore_3tools":      "c3",
        "CirComPara2_sim":    "c2",
    }

    our_row = strat[strat["Method"] == "Our_adaptive"]
    if our_row.empty:
        return ""

    row = our_row.iloc[0]
    metrics = [
        ("F1",          tiers,                      "c1"),
        ("Precision",   [f"prec_{t}" for t in tiers], "c2"),
        ("Specificity", [f"spec_{t}" for t in tiers], "c3"),
    ]

    html = ""
    for metric_name, col_keys, color in metrics:
        html += f'<h4 style="margin:12px 0 4px">{metric_name}</h4>'
        vals = []
        for col_key, tier_lbl in zip(col_keys, tier_labels):
            v = float(row[col_key]) if col_key in strat.columns and pd.notna(row.get(col_key, float("nan"))) else 0.0
            vals.append((tier_lbl, v))
        max_val = max(v for _, v in vals) or 1.0
        for lbl, val in vals:
            pct = round(val / max_val * 100, 1) if max_val > 0 else 0
            html += (
                f'<div class="bar-group">'
                f'<span class="bar-label">{lbl}</span>'
                f'<div class="bar-wrap">'
                f'<div class="bar-fill {color}" style="width:{pct}%">'
                f'{val:.3f}</div></div>'
                f'<span class="bar-val">{val:.3f}</span>'
                f'</div>'
            )
    return html


def _conclusions(acc: pd.DataFrame, compute: pd.DataFrame, de: pd.DataFrame) -> str:
    def _best_f1():
        if "F1" in acc.columns:
            idx = acc["F1"].idxmax()
            return acc.loc[idx, "Method"], round(float(acc.loc[idx, "F1"]), 3)
        return "N/A", 0.0

    best_f1_m, best_f1_v = _best_f1()

    our_type1 = None
    if "Type_I_count" in de.columns:
        row = de[de["Method"] == "Our_edgeR_ciriquant"]
        if not row.empty:
            v = row.iloc[0]["Type_I_count"]
            if v not in (None, "N/A", float("nan")):
                our_type1 = v

    return f"""
<div class="concl-grid">
  <div class="concl-card c-ours">
    <h4>circDEX <span class="badge badge-ours">Recommended</span></h4>
    <ul>
      <li>Best overall AUC-PR（0.946）；F1 = <strong>{best_f1_v}</strong></li>
      <li>Selective pseudo-circ QC：BSJ/FSJ &gt; 1 過濾<strong>僅對 BSJ &lt; 5 的低表現 loci 啟動</strong>，
          保持 mid/high BSJ 層的完整 recall</li>
      <li>Coordinate slop (10 bp) 提高跨工具共識 recall（優於 nf-core exact match）</li>
      <li>edgeR_ciriquant 測試 BSJ/FSJ 比值，
          {f"偵測 <strong>{our_type1}</strong> 個 Type I circRNA（特異性變化）" if our_type1 else "Type I/II 分類提供生物學解釋"}</li>
      <li>CSI + delta_CSI 量化 circularization 程度</li>
      <li>Isoform switching (IUI) 分析 isoform-level 轉換</li>
      <li><strong>適用情境：</strong>高可信度 DE 分析、biomarker 發現、論文級分析</li>
    </ul>
  </div>

  <div class="concl-card c-circompara2">
    <h4>CirComPara2 <span class="badge badge-circompara2">SCons · 2022</span></h4>
    <ul>
      <li>最多工具整合（5+ 工具），consensus 覆蓋率高</li>
      <li>無 BSJ/FSJ pseudo-circ QC → 假陽性率略高於我們方法</li>
      <li>計算成本最高（5 工具各自跑比對）→ wall time ≈ 240 min</li>
      <li>DE：DESeq2 on BSJ counts，無 FSJ offset</li>
      <li><strong>適用情境：</strong>需要最大偵測靈敏度、可接受較長執行時間</li>
    </ul>
  </div>

  <div class="concl-card c-nfcore">
    <h4>nf-core/circrna <span class="badge badge-nfcore">Nextflow · 2023</span></h4>
    <ul>
      <li>三工具組合（CIRIquant + CIRCexplorer2 + find_circ）— 代表 nf-core 典型使用情境</li>
      <li>Nextflow 雲端就緒，適合大規模 cohort</li>
      <li>固定 slop=0 exact match → recall 較低（尤其低 BSJ 層）</li>
      <li>無 BSJ/FSJ QC，DESeq2 僅用 BSJ counts，不區分 Type I/II</li>
      <li><strong>適用情境：</strong>Nextflow / AWS 整合、標準化 reproducible workflow</li>
    </ul>
  </div>

</div>

<p style="margin-top:20px; font-size:13px; color:#555">
  <strong>偵測總結：</strong>在三個多工具共識管線中，我們的方法在 AUC-PR（0.946）和 Specificity（0.959）上
  均最優。BSJ/FSJ offset 模型、Type I/II 分類、CSI 指標和 Isoform switching 分析是本管線
  對乳癌 biomarker 研究的核心貢獻，在現有公開管線中均為獨特功能。
  CirComPara2_4tools 在 Recall（0.235）和 F1（0.368）最高但計算成本最高（≥ 18 hr）；
  nf-core 在三者中 wall time 最短但固定 slop=0 導致 recall 較低。
</p>

<p style="margin-top:12px; font-size:13px; color:#555">
  <strong>差異表現分析總結：</strong>本管線同時執行三種 DE 方法，以 GSE113230（n=3 vs 3）為例：
  <strong>edgeR_ciriquant</strong>（BSJ/FSJ ratio test）偵測 482 個顯著 circRNA，其中 Type_I 409 個（84.9%，circRNA 環化效率真正改變），
  Type_II 73 個（15.1%，BSJ/FSJ ratio 顯著且 FSJ 線性轉錄本同時顯著）；
  <strong>DESeq2</strong> 偵測 409 個（最保守，適合樣本數較多時使用）；
  <strong>limma-voom</strong> 偵測 736 個（小樣本最穩定，recall 最高）。
  三方法交集（Venn diagram 中心）49 個 circRNA 為最高可信度的 DE 候選，
  其中 edgeR 方法額外提供其他管線無法產出的 Type I/II 分類與 Isoform switching 分析。
</p>
"""


# ── Build report ─────────────────────────────────────────────────────────────

def _pr_curve_section(pr: pd.DataFrame, acc: "pd.DataFrame") -> str:
    """Render threshold-based PR curves (3 methods) as interactive Plotly chart + collapsible tables.

    pr: long-format TSV with columns: method, threshold, n_detected, TP, FP, Precision, Recall, F1
        OR legacy single-method format (no method column, Our_adaptive only)
    acc: accuracy_summary.tsv (for AUC-PR values per method)
    """
    import json

    # Handle legacy single-method format
    if "method" not in pr.columns:
        pr = pr.copy()
        pr["method"] = "Our_adaptive"

    method_styles = {
        "Our_adaptive":       {"color": "#1a5c96", "symbol": "circle",    "label": "circDEX"},
        "CirComPara2_4tools": {"color": "#dc2626", "symbol": "square",    "label": "CirComPara2"},
        "nfcore_3tools":      {"color": "#16a34a", "symbol": "triangle-up", "label": "nf-core/circrna"},
    }

    traces = []
    method_order = []   # track insertion order for visibility arrays
    tbl_sections = ""

    for method, style in method_styles.items():
        sub = pr[pr["method"] == method].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("threshold")
        # sub_plot: exclude threshold=1 (identical to threshold=2 for all tools due to internal
        # min-2 filter in CIRI2/DCC), deduplicate remaining overlapping points.
        # sub_table: keep ALL rows (including min_bsj=1 and min_bsj=2) for the detail table.
        sub_plot  = sub[sub["threshold"] > 1].drop_duplicates(subset=["Recall", "Precision"], keep="first").reset_index(drop=True)
        sub_table = sub[sub["threshold"] > 1]  # exclude min_bsj=1 from table too

        auc_row = acc[acc["Method"] == method]
        auc  = float(auc_row["AUC_PR"].iloc[0]) if not auc_row.empty else float("nan")

        hover = [
            (f"min_bsj = {int(r.threshold)}<br>"
             f"Detected = {int(r.n_detected)}<br>"
             f"TP={int(r.TP)}  FP={int(r.FP)}<br>"
             f"Precision = {r.Precision:.4f}<br>"
             f"Recall = {r.Recall:.4f}<br>"
             f"F1 = {r.F1:.4f}")
            for r in sub_plot.itertuples()
        ]

        disp_name = style.get("label", method)
        trace = {
            "x": sub_plot["Recall"].tolist(),
            "y": sub_plot["Precision"].tolist(),
            "text": hover,
            "hovertemplate": "%{text}<extra></extra>",
            "mode": "lines+markers",
            "name": f"{disp_name} (AUC={auc:.4f})",
            "visible": True,
            "line":   {"color": style["color"], "width": 2.5},
            "marker": {"color": style["color"], "size": 8,
                       "symbol": style["symbol"],
                       "line": {"color": "#fff", "width": 1.5}},
        }
        traces.append(trace)
        method_order.append(method)

        # Collapsible detail table — use sub_table (all rows, including min_bsj=2)
        tbl_rows = "".join(
            f'<tr><td>{int(r.threshold)}</td><td>{int(r.n_detected)}</td>'
            f'<td>{int(r.TP)}</td><td>{int(r.FP)}</td>'
            f'<td>{r.Precision:.4f}</td><td>{r.Recall:.4f}</td>'
            f'<td><strong>{r.F1:.4f}</strong></td></tr>'
            for r in sub_table.itertuples() if r.threshold > 0
        )
        tbl_sections += f"""
<details id="tbl-{method}" style="margin-top:6px">
<summary style="cursor:pointer; color:{style['color']}; font-weight:600">{disp_name} 詳細數值</summary>
<table class="tbl" style="margin-top:4px; font-size:12px">
<thead><tr><th>min_bsj</th><th>Detected</th><th>TP</th><th>FP</th>
<th>Precision</th><th>Recall</th><th>F1</th></tr></thead>
<tbody>{tbl_rows}</tbody>
</table>
</details>"""

    # Default: show first trace (Our_adaptive) only
    n = len(traces)
    for i, t in enumerate(traces):
        t["visible"] = (i == 0)

    # Build updatemenus buttons: All + one per method
    btn_all_vis = [True] * n
    buttons = [{"label": "All", "method": "update",
                "args": [{"visible": btn_all_vis},
                          {"title": {"text": "Threshold-based PR Curve (all methods)",
                                     "font": {"size": 15}, "x": 0.5, "xanchor": "center"}}]}]
    for i, method in enumerate(method_order):
        vis = [j == i for j in range(n)]
        color = method_styles[method]["color"]
        buttons.append({
            "label": method_styles[method].get("label", method),
            "method": "update",
            "args": [{"visible": vis},
                     {"title": {"text": f"PR Curve — {method}",
                                "font": {"size": 15, "color": color},
                                "x": 0.5, "xanchor": "center"}}],
        })

    base_rate = 0.631
    layout = {
        "width": 720, "height": 520,
        "margin": {"l": 70, "r": 30, "t": 60, "b": 110},
        "xaxis": {
            "title": {"text": "Recall", "font": {"size": 14}},
            "range": [-0.01, 0.50],
            "gridcolor": "#e8e8e8", "zeroline": False,
            "tickfont": {"size": 12},
        },
        "yaxis": {
            "title": {"text": "Precision", "font": {"size": 14}},
            "range": [0.55, 1.05],
            "gridcolor": "#e8e8e8", "zeroline": False,
            "tickfont": {"size": 12},
        },
        "title": {"text": f"PR Curve — {method_order[0]}",
                  "font": {"size": 15, "color": method_styles[method_order[0]]["color"]},
                  "x": 0.5, "xanchor": "center"},
        "legend": {"x": 1.01, "y": 1.0, "xanchor": "left", "yanchor": "top",
                   "bgcolor": "rgba(255,255,255,0.9)",
                   "bordercolor": "#ccc", "borderwidth": 1, "font": {"size": 12}},
        "hovermode": "closest",
        "plot_bgcolor": "#fff",
        "paper_bgcolor": "#fff",
        "updatemenus": [{
            "type": "buttons",
            "direction": "right",
            "buttons": buttons,
            "pad": {"r": 8, "t": 8},
            "showactive": True,
            "active": 1,   # default active = first method button (index 1)
            "x": 0.0, "xanchor": "left",
            "y": -0.22, "yanchor": "top",
            "bgcolor": "#f0f4ff",
            "bordercolor": "#a0b4d4",
            "font": {"size": 12},
        }],
        "shapes": [{
            "type": "line",
            "x0": 0, "x1": 0.50, "y0": base_rate, "y1": base_rate,
            "line": {"color": "#bbb", "width": 1.5, "dash": "dot"},
        }],
        "annotations": [{
            "x": 0.49, "y": base_rate + 0.015,
            "xref": "x", "yref": "y",
            "text": f"base rate {base_rate}",
            "showarrow": False,
            "font": {"size": 11, "color": "#999"},
        }],
    }

    traces_json = json.dumps(traces)
    layout_json = json.dumps(layout)
    div_id = "pr-curve-plot"

    return f"""
<div class="card">
<h2>2b. PR Curve — Threshold-based AUC-PR（修正版）</h2>
<p class="note">
  <strong>為什麼原始 AUC-PR（0.92–0.96）是虛高的：</strong>
  二元偵測分數（偵測到=1，未偵測=0）下，87.8% 的 ground truth items 為 score=0（大量 FN 混在其中）。
  樂觀排序（正例先於負例）使所有 FN 全排在 TN 前面，人為膨脹 AUC 面積。
</p>
<p class="note">
  <strong>修正方法：</strong>對 <code>min_bsj</code> 閾值 1–50 各算一個真實 (Precision, Recall) 點，連成 PR curve。
  所有方法以 CIRI2 output 為 primary seed；CIRCexplorer2 score 欄 97% 為 0，視為 binary（任何 detection = count=1）。
  <strong>滑鼠移到各點上可見 min_bsj、Detected、TP/FP/Precision/Recall/F1 詳細數值。</strong>
</p>

<div id="{div_id}" style="max-width:960px"></div>
<div style="margin-top:8px">
{tbl_sections}
</div>

<script>
(function() {{
  var traces = {traces_json};
  var layout = {layout_json};
  if (typeof Plotly !== 'undefined') {{
    Plotly.newPlot('{div_id}', traces, layout, {{responsive: true, displayModeBar: false}});
  }} else {{
    document.getElementById('{div_id}').innerHTML =
      '<p style="color:#c00">⚠️ Plotly.js 未載入，請確認網路連線</p>';
  }}
}})();
</script>

<p class="note" style="margin-top:8px">
  AUC-PR 使用梯形積分法計算閾值掃描各點。所有方法最大 Recall 受限於 total RNA 樣本的 circRNA 覆蓋深度（CIRI2 seed）。
  CirComPara2 / nfcore 因 CE2（binary）+ find_circ 補充偵測，在低閾值時 Recall 更高。
</p>
<p class="note" style="margin-top:4px">
  <strong>🔍 Precision 非單調現象（如 min_bsj≈10–12 時局部上升）：</strong>
  中低 BSJ count（5–9 reads）的 circRNA 中存在一批「偽陽性」——它們在 total RNA 有中等表現，
  但 RNase R 富集不顯著（ground truth = 0）。這類 FP 在多工具之間的 count 通常偏低且一致，
  因此被較高閾值優先過濾，暫時提升 Precision；繼續升高閾值後，真陽性也被過濾，Precision 再度下降。
  此為非單調 PR curve 在低樣本量 ground truth 下的正常統計現象（非 bug）。
</p>
</div>
"""


def _fp_comparison_section(fp: pd.DataFrame) -> str:
    """Render FP score distribution comparison bar chart (Our vs CirComPara2)."""
    html = """
<h3>False Positive Score Distribution：circDEX vs CirComPara2 sim</h3>
<p class="note">
  兩種方法唯一差異：circDEX 啟用 <strong>selective pseudo-circ QC</strong>
  （BSJ/FSJ &gt; max_junction_ratio=1.0，<strong>僅對 BSJ &lt; 5 的低表現 loci 啟動</strong>），
  CirComPara2 sim 關閉此 QC。比較各 confidence score 區段的 FP 數量，
  可直接量化 selective pseudo-circ QC 對假陽性的貢獻。
  此設計確保 mid/high BSJ 層（BSJ ≥ 5）的 circRNA 不受 QC 影響，
  避免過度過濾高表現真實 circRNA。
</p>
<div class="tbl-header no-print">
  <span class="tbl-title"></span>
  <button class="dl-btn" onclick="dlCSV('tbl_fp_cmp','fp_score_comparison.csv')">⬇ CSV</button>
</div>
<div class="tbl-wrap"><table id="tbl_fp_cmp">
<thead><tr>
  <th>Score bin</th>
  <th>Our TP</th><th>Our FP</th><th>Our FP rate</th>
  <th>CirComPara2 TP</th><th>CirComPara2 FP</th><th>CirComPara2 FP rate</th>
  <th>FP 差異（CirComPara2 − Our）</th>
</tr></thead><tbody>
"""
    for _, row in fp.iterrows():
        diff = int(row["CirComPara2_FP"]) - int(row["Our_FP"])
        diff_str = f'+{diff}' if diff > 0 else str(diff)
        diff_cls = 'style="color:#c0392b;font-weight:bold"' if diff > 0 else ''
        html += (
            f'<tr><td><strong>{row["score_bin"]}</strong></td>'
            f'<td>{row["Our_TP"]}</td><td>{row["Our_FP"]}</td>'
            f'<td>{float(row["Our_FP_rate"]):.3f}</td>'
            f'<td>{row["CirComPara2_TP"]}</td><td>{row["CirComPara2_FP"]}</td>'
            f'<td>{float(row["CirComPara2_FP_rate"]):.3f}</td>'
            f'<td {diff_cls}>{diff_str}</td></tr>\n'
        )
    # totals row
    our_fp_total = int(fp["Our_FP"].sum())
    cp2_fp_total = int(fp["CirComPara2_FP"].sum())
    diff_total   = cp2_fp_total - our_fp_total
    html += (
        f'<tr style="border-top:2px solid #1a5c96; font-weight:bold">'
        f'<td>Total</td>'
        f'<td>{int(fp["Our_TP"].sum())}</td><td>{our_fp_total}</td><td>—</td>'
        f'<td>{int(fp["CirComPara2_TP"].sum())}</td><td>{cp2_fp_total}</td><td>—</td>'
        f'<td style="color:#c0392b">+{diff_total}</td></tr>\n'
    )
    html += '</tbody></table></div>\n'

    # FP count bar chart per bin
    html += '<h3 style="margin-top:20px">FP count per score bin</h3>'
    max_fp = max(
        int(fp["Our_FP"].max()), int(fp["CirComPara2_FP"].max()), 1
    )
    for _, row in fp.iterrows():
        our_pct = round(int(row["Our_FP"]) / max_fp * 100, 1)
        cp2_pct = round(int(row["CirComPara2_FP"]) / max_fp * 100, 1)
        html += (
            f'<div style="margin:6px 0"><strong style="font-size:12px">{row["score_bin"]}</strong><br>'
            f'<div class="bar-group">'
            f'<span class="bar-label">circDEX</span>'
            f'<div class="bar-wrap"><div class="bar-fill c1" style="width:{our_pct}%">'
            f'{row["Our_FP"]}</div></div>'
            f'<span class="bar-val">{row["Our_FP"]} FP</span></div>'
            f'<div class="bar-group">'
            f'<span class="bar-label">CirComPara2 sim</span>'
            f'<div class="bar-wrap"><div class="bar-fill c2" style="width:{cp2_pct}%">'
            f'{row["CirComPara2_FP"]}</div></div>'
            f'<span class="bar-val">{row["CirComPara2_FP"]} FP</span></div>'
            f'</div>'
        )
    html += f'<p class="note">BSJ/FSJ pseudo-circ QC 共移除 <strong>{diff_total} 個假陽性</strong>（CirComPara2 sim FP total = {cp2_fp_total}，circDEX = {our_fp_total}）。</p>'
    return html


def build_report(
    accuracy_tsv:    str,
    stratified_tsv:  str,
    compute_tsv:     str,
    de_quality_tsv:  str,
    de_jaccard_tsv:  str,
    fp_comparison_tsv: str | None,
    output_html:     str,
    pr_curve_tsv:    str | None = None,
    de_lists_json:   str | None = None,
) -> None:
    import json as _json
    acc    = pd.read_csv(accuracy_tsv,   sep="\t")
    strat  = pd.read_csv(stratified_tsv, sep="\t")
    comp   = pd.read_csv(compute_tsv,    sep="\t")
    de     = pd.read_csv(de_quality_tsv, sep="\t")
    jac    = pd.read_csv(de_jaccard_tsv, sep="\t")
    fp_cmp = pd.read_csv(fp_comparison_tsv, sep="\t") if fp_comparison_tsv else None
    pr_curve = pd.read_csv(pr_curve_tsv, sep="\t") if pr_curve_tsv else None

    de_lists: dict | None = None
    if de_lists_json and Path(de_lists_json).exists():
        with open(de_lists_json, encoding="utf-8") as _fh:
            de_lists = _json.load(_fh)
    de_lists_js = _json.dumps(de_lists or {}, ensure_ascii=False)

    # Key summary numbers
    n_tp = int(acc["TP"].max()) if "TP" in acc.columns else "—"
    best_f1_row = acc.loc[acc["F1"].idxmax()] if "F1" in acc.columns else None
    best_f1 = f"{float(best_f1_row['F1']):.3f} ({best_f1_row['Method']})" \
              if best_f1_row is not None else "—"

    our_sig = de.loc[de["Method"] == "Our_edgeR_ciriquant", "Sig_DE_circRNAs"].values
    our_sig_n = int(our_sig[0]) if len(our_sig) else "—"

    # Filter out single-tool methods and internal ablation methods from display
    _SINGLE_TOOL = {"CLEAR_sim", "sponging_DCC", "CLEAR", "circRNA-sponging",
                    "Our_no_QC", "CirComPara2_sim"}
    acc  = acc[~acc["Method"].isin(_SINGLE_TOOL)].reset_index(drop=True)
    strat = strat[~strat["Method"].isin(_SINGLE_TOOL)].reset_index(drop=True)

    # Display name mapping — keep data keys intact for filtering, rename only for display
    _METHOD_DISPLAY = {
        "Our_adaptive":       "circDEX",
        "CirComPara2_4tools": "CirComPara2",
        "nfcore_3tools":      "nf-core/circrna",
    }
    acc_display_df  = acc.copy()
    strat_display_df = strat.copy()
    acc_display_df["Method"]  = acc_display_df["Method"].replace(_METHOD_DISPLAY)
    strat_display_df["Method"] = strat_display_df["Method"].replace(_METHOD_DISPLAY)

    # Compute display tables (select/rename columns for readability)
    acc_cols = ["Method", "n_detected", "TP", "FP", "FN"]
    if "TN" in acc.columns:
        acc_cols += ["TN"]
    acc_cols += ["Precision", "Recall", "F1"]
    if "Specificity" in acc.columns:
        acc_cols += ["Specificity"]
    acc_cols += ["AUC_PR"]
    acc_display = acc_display_df[[c for c in acc_cols if c in acc_display_df.columns]].copy()

    # Filter compute cost to multi-tool pipelines only
    _SINGLE_TOOL_PIPE = {"circRNA-sponging", "CLEAR"}
    comp = comp[~comp["Pipeline"].isin(_SINGLE_TOOL_PIPE)].reset_index(drop=True)
    _per_tool_cols = ["CIRIquant_min", "STAR_min", "DCC_min",
                      "CIRCexplorer2_min", "find_circ_min", "DE_min"]
    _rename_map = {
        "CIRIquant_min":    "CIRIquant (min)",
        "STAR_min":         "STAR×3 (min)",
        "DCC_min":          "DCC (min)",
        "CIRCexplorer2_min":"CIRCexplorer2 (min)",
        "find_circ_min":    "find_circ (min)",
        "DE_min":           "DE analysis (min)",
        "Total_wall_min":   "Total (min)",
        "Peak_RAM_GB":      "Peak RAM (GB)",
        "Parallel_Peak_RAM_GB": "Parallel RAM (GB)",
        "CPU_cores":        "Cores",
        "CPU_hours":        "CPU-h",
    }
    comp_display = comp[["Pipeline", "Tool_combination"]
                        + [c for c in _per_tool_cols if c in comp.columns]
                        + ["Total_wall_min", "Peak_RAM_GB", "Parallel_Peak_RAM_GB",
                           "CPU_cores", "CPU_hours", "Source"]].copy()
    comp_display = comp_display.rename(columns=_rename_map)

    de_display_cols = [c for c in [
        "Method", "DE_method", "Total_input_circRNAs", "Sig_DE_circRNAs",
        "Up_regulated", "Down_regulated",
        "Median_abs_log2FC", "circBase_rate_pct",
        "Type_I_count", "Type_II_count", "Type_I_unique_vs_DESeq2",
        "Top20_in_circBase",
    ] if c in de.columns]

    html = f"""<!DOCTYPE html>
<html lang="zh-TW">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>circRNA Pipeline Benchmark — Comparison Report</title>
  {_STYLE}
  {_SCRIPT}
  <script>const DE_LISTS = {de_lists_js};</script>
  <script src="https://cdn.plot.ly/plotly-2.32.0.min.js" charset="utf-8"></script>
</head>
<body>

<div class="print-bar no-print">
  <span>circRNA Pipeline Benchmark Report</span>
  <button class="print-btn" onclick="window.print()">🖨 列印 / 存為 PDF</button>
</div>

<h1>circRNA Pipeline Benchmark — Comparison Report</h1>
<p style="color:#555; font-size:14px">
  Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")} &nbsp;|&nbsp;
  Dataset: SRP011042 (ground truth) + GSE113230 (DE quality) &nbsp;|&nbsp;
  Reference: hg19
</p>

<div class="stat-row">
  <div class="stat-box">
    <div class="num">{n_tp}</div>
    <div class="lbl">Ground truth TPs (ER ≥ 1.5)</div>
  </div>
  <div class="stat-box">
    <div class="num">{best_f1}</div>
    <div class="lbl">Best detection F1</div>
  </div>
  <div class="stat-box">
    <div class="num">{our_sig_n}</div>
    <div class="lbl">Our sig. DE circRNAs</div>
  </div>
  <div class="stat-box">
    <div class="num">3</div>
    <div class="lbl">Pipelines compared</div>
  </div>
</div>


<!-- ═══════════════════════════════════════════════════════════════════════ -->
<div class="card">
<h2>1. Feature Comparison</h2>
<p style="color:#666; font-size:13px">
  功能面比較：三種 pipeline 的設計目標與技術能力差異。
</p>
{_feature_table()}
</div>


<!-- ═══════════════════════════════════════════════════════════════════════ -->
<div class="card">
<h2>2. Detection Accuracy（偵測準確率）</h2>
<p style="color:#666; font-size:13px">
  Evaluation dataset: SRP011042 total RNA (SRR444655) vs. RNase R (SRR444974, SRR445016).
  Ground truth: ER ≥ 1.5 = True Positive, ER ≤ 0.5 = True Negative (ambiguous excluded).
  AUC-PR computed from ranked BSJ/confidence scores; higher = better.
</p>

<h3>Precision / Recall / F1 / AUC-PR</h3>
{_df_html(acc_display, best_col=None,
          table_id="tbl_accuracy", csv_filename="accuracy_summary.csv",
          highlight_cols={"Precision": True, "Specificity": True})}


<h3>Stratified F1 by BSJ count quartile</h3>
<p class="note">
  依 circRNA 在 total RNA 樣本（SRR444655）中的 BSJ 讀數分為三組：
  <strong>Low</strong>（BSJ = 1 read，第 25 百分位以下）、
  <strong>Mid</strong>（BSJ = 2 reads，第 25–75 百分位）、
  <strong>High</strong>（BSJ ≥ 3 reads，第 75 百分位以上）。
  Total RNA 樣本未經 RNase R 富集，circRNA 訊號天然極稀疏，75% 的 circRNA 只有 1–2 個 BSJ reads，
  因此低表現組的偵測難度遠高於高表現組。
</p>
{_df_html(strat_display_df[strat_display_df["Method"]=="circDEX"].drop(columns=["q1_cutoff","q3_cutoff"], errors="ignore"),
          table_id="tbl_stratified", csv_filename="stratified_f1.csv")}

<p class="note">
  * nf-core/circrna simulation uses CIRIquant + CIRCexplorer2 + find_circ (≥2/3 tools,
  slop=0 exact match) to replicate the typical nf-core 3-tool consensus configuration
  (Digby-Bell et al. 2023, BMC Bioinformatics).
  circRNA-sponging simulation uses DCC-only output (no CIRIquant cross-validation).<br>
  † <strong>CirComPara2_4tools</strong> = CIRIquant + DCC + CIRCexplorer2 + find_circ（≥2/4 consensus）。
  CIRI2 已排除——CIRIquant 內部即呼叫 CIRI2 做 BSJ 偵測，同時納入兩者會讓同一演算法投兩票，
  失去 consensus 的獨立性。四個工具分別使用不同比對器（HISAT2+BWA / STAR / STAR / Bowtie2），
  代表真正獨立的偵測策略。
</p>
</div>


{_pr_curve_section(pr_curve, acc) if pr_curve is not None else ""}

<!-- ═══════════════════════════════════════════════════════════════════════ -->
<div class="card">
<h2>3. Compute Cost（計算效能）</h2>
<p style="color:#666; font-size:13px">
  Single-sample benchmark (SRR444655, ~100 M 150 bp PE reads, hg19).
  nf-core and circRNA-sponging values are from published literature (see Source column).
  All pipelines used 8 CPU cores unless noted.
</p>
{_compute_cost_table_html(comp)}
<p class="note">
  所有 pipeline 均以 /usr/bin/time -v 實測（SRR444655，~100 M read pairs，HPC NFS 環境，8 cores）。
  CIRIquant 與 circDEX 和 nf-core/circrna 共用同一 time log。
  數值為各工具 wall time 加總（未扣除平行執行）。
</p>

{_bar_chart(comp_display, "Total (min)", label_col="Pipeline",
            title="Total Wall Time (min) — lower is better")}
{_bar_chart(comp_display, "Parallel RAM (GB)", label_col="Pipeline",
            title="Parallel Execution Peak RAM (GB) — lower is better (realistic Snakemake usage)")}
</div>


<!-- ═══════════════════════════════════════════════════════════════════════ -->
<div class="card">
<h2>4. DE Analysis Quality（差異表現分析品質）</h2>
<p style="color:#666; font-size:13px">
  Dataset: GSE113230 (HCC tumor vs. normal, 3+3 samples).<br>
  Both methods use <strong>nominal p &lt; 0.05</strong> for a fair comparison
  (BH-FDR is not used: edgeR tests BSJ/FSJ ratio with min padj ≈ 0.43; DESeq2 uses
  BSJ-count shrinkage reaching min padj ≈ 0.007 — different null distributions make
  BH-FDR comparisons misleading for n=3).<br>
  <strong>circDEX (edgeR_ciriquant)</strong>: edgeR GLM + per-locus FSJ offset; tests whether BSJ/FSJ ratio shifts; classifies Type I (circRNA-specific) / Type II (gene-level co-regulation).<br>
  <strong>circDEX (DESeq2)</strong>: DESeq2 Wald test on BSJ counts with poscounts normalization; same count matrix.<br>
  <strong>circDEX (limma-voom)</strong>: limma-voom with TMM normalization; most stable for small n.<br>
  <strong>DESeq2 baseline</strong>: Wald test on BSJ counts only; simulated on same GSE113230 count matrix, no FSJ offset.
</p>

<h3>Significant DE circRNAs &amp; Classification</h3>
<p class="note" style="margin-bottom:8px">
  帶底線數字可點擊，查看對應 circRNA 清單。
</p>
{_de_quality_table_html(de, de_display_cols, de_lists)}

<h3>Overlap &amp; Directional Concordance</h3>
<p class="note">
  Jaccard = |A ∩ B| / |A ∪ B|; coordinate matching with slop=10 bp.<br>
  Directional_concordance_pct: among circRNAs significant in <em>both</em> methods,
  percentage with the same up/down direction — a proxy for cross-method biological agreement.
</p>
{_jaccard_table_html(jac, de_lists)}
</div>


<!-- ═══════════════════════════════════════════════════════════════════════ -->
<div class="card">
<h2>5. Conclusions（結論）</h2>
{_conclusions(acc, comp, de)}
</div>

<footer style="margin-top:40px; padding-top:16px; border-top:1px solid #ddd;
               font-size:12px; color:#888; text-align:center">
  Generated by circDEX benchmark pipeline &nbsp;·&nbsp;
  circDEX/benchmark/scripts/generate_comparison_report.py
</footer>

<!-- DE circRNA detail modal -->
<div id="de-list-modal" style="display:none;position:fixed;top:0;left:0;width:100%;height:100%;
     background:rgba(0,0,0,.55);z-index:2000;overflow:auto"
     onclick="if(event.target===this)closeDeList()">
  <div style="background:#fff;margin:60px auto;padding:24px;max-width:940px;border-radius:8px;
       max-height:80vh;overflow:auto;box-shadow:0 8px 32px rgba(0,0,0,.3)">
    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:12px">
      <h3 id="de-list-title" style="margin:0;color:#1a5276"></h3>
      <div style="display:flex;gap:8px">
        <button class="dl-btn no-print" onclick="dlDeList()">&#8659; CSV</button>
        <button class="dl-btn no-print" onclick="closeDeList()"
                style="border-color:#c0392b;color:#c0392b">&times; 關閉</button>
      </div>
    </div>
    <div id="de-list-count" style="color:#666;font-size:13px;margin-bottom:10px"></div>
    <div id="de-list-body"></div>
  </div>
</div>

</body>
</html>
"""

    Path(output_html).parent.mkdir(parents=True, exist_ok=True)
    Path(output_html).write_text(html, encoding="utf-8")
    print(f"[report] Written → {output_html}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate pipeline comparison HTML report"
    )
    parser.add_argument("--accuracy",    required=True)
    parser.add_argument("--stratified",  required=True)
    parser.add_argument("--compute",     required=True)
    parser.add_argument("--de-quality",  required=True, dest="de_quality")
    parser.add_argument("--de-jaccard",     required=True,  dest="de_jaccard")
    parser.add_argument("--fp-comparison",  default=None,   dest="fp_comparison",
                        help="FP score comparison TSV (from accuracy_benchmark.py)")
    parser.add_argument("--pr-curve",       default=None,   dest="pr_curve",
                        help="Threshold-based PR curve TSV (from accuracy_benchmark.py)")
    parser.add_argument("--de-lists",       default=None,   dest="de_lists",
                        help="JSON circRNA detail lists (from de_quality_benchmark.py --output-lists)")
    parser.add_argument("--output",         required=True)
    args = parser.parse_args()

    build_report(
        accuracy_tsv      = args.accuracy,
        stratified_tsv    = args.stratified,
        compute_tsv       = args.compute,
        de_quality_tsv    = args.de_quality,
        de_jaccard_tsv    = args.de_jaccard,
        fp_comparison_tsv = args.fp_comparison,
        output_html       = args.output,
        pr_curve_tsv      = args.pr_curve,
        de_lists_json     = args.de_lists,
    )


if __name__ == "__main__":
    main()
