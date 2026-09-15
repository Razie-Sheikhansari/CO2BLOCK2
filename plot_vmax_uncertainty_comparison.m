% Interactive HTML chart (log y-axis): dynamic capacity lines per well count
% and static compressibility-based ceiling lines per mode.
% Static ceilings are horizontal dotted lines, same color as matching dynamic line.
%
% Usage:
%   plot_static_capacity_comparison()           % loads CO2BLOCK_results.mat
%   plot_static_capacity_comparison(OUT, R)
%   plot_static_capacity_comparison(OUT, R, 'out.html')

function plot_static_capacity_comparison(OUT, R, html_out)

%% 0. Load
if nargin < 1 || isempty(OUT)
    loaded = load('CO2BLOCK_results.mat', 'OUT', 'R');
    OUT    = loaded.OUT;
    R      = loaded.R;
    fprintf('Loaded results from CO2BLOCK_results.mat\n');
elseif nargin < 2 || isempty(R)
    loaded = load('CO2BLOCK_results.mat', 'R');
    R      = loaded.R;
end
if nargin < 3 || isempty(html_out)
    html_out = 'plot_static_capacity.html';
end

%% 1. Flags
is_fault      = (OUT.center_point_type == "Fault");
stress_active = isfield(OUT,'stress_uncertain')  && OUT.stress_uncertain  == "yes";
fault_active  = isfield(OUT,'fault_uncertainty') && OUT.fault_uncertainty == "yes";
perm_active   = isfield(OUT,'perm_uncertain')    && logical(OUT.perm_uncertain);
any_uncertain = stress_active || fault_active || perm_active;
is_det        = ~any_uncertain;
% perm only: deltap_cr doesn't vary, so one static ceiling applies to all modes
perm_only_static = false;

well_list   = OUT.well_list(:)';
d_list      = OUT.d_list(:)';
x_grid_list = OUT.x_grid_list(:)';
y_grid_list = OUT.y_grid_list(:)';
nw          = numel(well_list);
nd          = numel(d_list);

%% 2. Static capacity
%  M_static [Gt] = V_bulk [m3] * compr [Pa-1] * DeltaP_cr [Pa] * dens_c [kg/m3] / 1e12
V_bulk = R.area_res * 1e6 * R.thick;          % km2 -> m2, then x thickness [m]
gt     = @(dc_MPa) V_bulk * R.compr * (dc_MPa * 1e6) * R.dens_c / 1e12;  % [Gt]

if is_det
    % deterministic
    dc_det_arr = OUT.deltap_cr_det;
    dc_fault   = min(dc_det_arr(dc_det_arr > 0), [], 'omitnan');  % most critical fault [MPa]
    dc_ref     = safe_scalar(OUT, 'deltap_cr_ref', dc_fault);     % central point

    if is_fault
        % two ceilings: reservoir ref point and fault constraint
        static_dc   = [dc_ref,       dc_fault      ];
        static_Gt   = [gt(dc_ref),   gt(dc_fault)  ];
        static_lbls = {'Static Det', 'Static Fault Det'};
        static_cols = {'#444455',    '#7a7a8a'       };
    else
        % single ceiling
        static_dc   = dc_fault;
        static_Gt   = gt(dc_fault);
        static_lbls = {'Static Det'};
        static_cols = {'#444455'  };
    end

else
    % MC: any uncertainty active
    samps = OUT.deltap_cr_samples;
    if iscell(samps)
        % fault mode: each cell = one fault's sample array; take min across faults
        samp_mat    = cell2mat(samps(:)');
        dc_per_real = min(samp_mat, [], 2, 'omitnan');
    elseif ~isempty(samps)
        dc_per_real = samps(:);
    else
        dc_per_real = [];
    end
    if ~isempty(dc_per_real)
        dc_per_real(dc_per_real <= 0) = NaN;
    end

    % deterministic reference
    dc_det_arr = OUT.deltap_cr_det;
    dc_det     = min(dc_det_arr(dc_det_arr > 0), [], 'omitnan');

    % if only perm_uncertain=yes, deltap_cr_samples is empty (perm doesn't affect deltap_cr)
    have_dc_mc = ~isempty(dc_per_real) && any(isfinite(dc_per_real));

    if have_dc_mc
        % stress or fault uncertain: deltap_cr varies across realizations
        pcts = prctile(dc_per_real, [90, 50, 10]);  % [P90, P50, P10]

        if is_fault
            % 5 ceilings: Det, Fault Det, P90, P50, P10
            dc_ref = safe_scalar(OUT, 'deltap_cr_ref', dc_det);
            static_dc   = [dc_ref,       dc_det,              pcts(1),       pcts(2),       pcts(3)     ];
            static_Gt   = arrayfun(gt, static_dc);
            static_lbls = {'Static Det', 'Static Fault Det',  'Static P90',  'Static P50',  'Static P10'};
            static_cols = {'#444455',    '#7a7a8a',            '#1a8c38',     '#e09010',     '#d43030'   };
        else
            % 4 ceilings: Det, P90, P50, P10
            static_dc   = [dc_det,       pcts(1),       pcts(2),       pcts(3)     ];
            static_Gt   = arrayfun(gt, static_dc);
            static_lbls = {'Static Det', 'Static P90',  'Static P50',  'Static P10'};
            static_cols = {'#444455',    '#1a8c38',     '#e09010',     '#d43030'   };
        end

    else
        % perm only: one shared ceiling for all modes
        static_dc   = dc_det;
        static_Gt   = gt(dc_det);
        static_lbls = {'Static bound - all modes'};
        static_cols = {'#6644aa'};
        perm_only_static = true;
    end
end

n_static = numel(static_Gt);
fprintf('Static capacities:\n');
for si = 1:n_static
    fprintf('  %-18s  DeltaP_cr = %.3f MPa  ->  %.4f Gt\n', ...
        static_lbls{si}, static_dc(si), static_Gt(si));
end

%% 3. Dynamic capacity lines
% colors match corresponding static ceiling (tooltip uses color to pair them)
dyn_lbl  = {};
dyn_col  = {};
dyn_DATA = {};
dyn_dash = {};

% Det
det_data = safe_field(OUT, 'V_M_ref', safe_field(OUT, 'V_M_det', []));
if ~isempty(det_data) && any(isfinite(det_data(:)))
    dyn_lbl{end+1}  = 'Det dynamic';
    dyn_col{end+1}  = '#444455';
    dyn_DATA{end+1} = det_data;
    dyn_dash{end+1} = '[8,4]';
end

% Fault Det
if is_fault
    fd = safe_field(OUT, 'V_M_fault_det', []);
    if ~isempty(fd) && any(isfinite(fd(:)))
        dyn_lbl{end+1}  = 'Fault Det dynamic';
        dyn_col{end+1}  = '#7a7a8a';
        dyn_DATA{end+1} = fd;
        dyn_dash{end+1} = '[4,3]';
    end
end

% MC: P10, P50, P90
if any_uncertain
    jcols = {'#d43030','#e09010','#1a8c38'};
    if is_fault
        sets = {safe_field(OUT,'V_M_joint_fault_P10', safe_field(OUT,'V_M_joint_P10',[])), ...
                safe_field(OUT,'V_M_joint_fault_P50', safe_field(OUT,'V_M_joint_P50',[])), ...
                safe_field(OUT,'V_M_joint_fault_P90', safe_field(OUT,'V_M_joint_P90',[]))};
    else
        sets = {safe_field(OUT,'V_M_joint_P10',[]), ...
                safe_field(OUT,'V_M_joint_P50',[]), ...
                safe_field(OUT,'V_M_joint_P90',[])};
    end
    jlbls = {'P10 dynamic','P50 dynamic','P90 dynamic'};
    for k = 1:3
        if ~isempty(sets{k}) && any(isfinite(sets{k}(:)))
            dyn_lbl{end+1}  = jlbls{k};
            dyn_col{end+1}  = jcols{k};
            dyn_DATA{end+1} = sets{k};
            dyn_dash{end+1} = '[]';
        end
    end
end

n_dyn = numel(dyn_lbl);
if n_dyn == 0
    warning('No dynamic data found. Nothing to plot.'); return;
end

%% 4. Best scenario per well count
n_sc     = nw * nd;
scen_nw  = zeros(n_sc,1);
scen_nx  = zeros(n_sc,1);
scen_ny  = zeros(n_sc,1);
scen_dkm = zeros(n_sc,1);
dyn_mat  = NaN(n_sc, n_dyn);
sc = 0;
for wi = 1:nw
    for di = 1:nd
        sc = sc + 1;
        scen_nw(sc)  = well_list(wi);
        scen_nx(sc)  = x_grid_list(wi);
        scen_ny(sc)  = y_grid_list(wi);
        scen_dkm(sc) = d_list(di);
        for li = 1:n_dyn
            M = dyn_DATA{li};
            if isempty(M) || wi > size(M,1) || di > size(M,2), continue; end
            v = M(wi,di);
            if isfinite(v) && v >= 0, dyn_mat(sc,li) = v; end
        end
    end
end

valid    = any(isfinite(dyn_mat), 2);
dyn_mat  = dyn_mat(valid,:);
scen_nw  = scen_nw(valid);
scen_nx  = scen_nx(valid);
scen_ny  = scen_ny(valid);
scen_dkm = scen_dkm(valid);

% keep best scenario per unique well count (rank by Det line 1)
uniq_nw  = unique(scen_nw);
keep_idx = zeros(numel(uniq_nw), 1);
for k = 1:numel(uniq_nw)
    rows = find(scen_nw == uniq_nw(k));
    crit = dyn_mat(rows, 1);
    if all(isnan(crit)), crit = max(dyn_mat(rows,:), [], 2); end
    [~, best]   = max(crit);
    keep_idx(k) = rows(best);
end

dyn_mat  = dyn_mat(keep_idx, :);
scen_nw  = scen_nw(keep_idx);
scen_nx  = scen_nx(keep_idx);
scen_ny  = scen_ny(keep_idx);
scen_dkm = scen_dkm(keep_idx);
N        = numel(keep_idx);

%% 5. Chart.js dataset strings
ds_parts = {};

% dynamic datasets
for li = 1:n_dyn
    vals = dyn_mat(:, li);
    vstr = strjoin(arrayfun(@(v) ifelse(isfinite(v), sprintf('%.4f',v), 'null'), ...
        vals, 'UniformOutput', false), ',');
    is_d = contains(dyn_lbl{li}, 'Det');
    bw   = ifelse(is_d, '1.5', '2.5');
    pr   = ifelse(is_d, '2.5', '4.5');
    phr  = ifelse(is_d, '5',   '8');
    ds_parts{end+1} = sprintf([...
'{label:"%s",data:[%s],borderColor:"%s",'...
'backgroundColor:"%s33",borderWidth:%s,'...
'pointRadius:%s,pointHoverRadius:%s,'...
'pointBackgroundColor:"%s",pointBorderColor:"#fff",'...
'pointBorderWidth:1.5,borderDash:%s,tension:0.25,fill:false}'], ...
        dyn_lbl{li}, vstr, dyn_col{li}, dyn_col{li}, bw, pr, phr, ...
        dyn_col{li}, dyn_dash{li});
end

% static datasets (same value at every x point)
for si = 1:n_static
    sv   = static_Gt(si);
    vstr = strjoin(repmat({sprintf('%.4f', sv)}, 1, N), ',');
    ds_parts{end+1} = sprintf([...
'{label:"%s  (DeltaPcr=%.2f MPa, %.4f Gt)",'...
'data:[%s],borderColor:"%s",'...
'backgroundColor:"transparent",borderWidth:1.5,'...
'pointRadius:0,pointHoverRadius:0,'...
'borderDash:[3,3],tension:0,fill:false}'], ...
        static_lbls{si}, static_dc(si), sv, vstr, static_cols{si});
end

datasets_js = strjoin(ds_parts, ',');

% JS metadata arrays
xlabs       = strjoin(arrayfun(@(n) sprintf('"%d"',n), scen_nw, 'UniformOutput',false), ',');
wells_js    = ['[' strjoin(arrayfun(@num2str, scen_nw,  'UniformOutput',false), ',') ']'];
nx_js       = ['[' strjoin(arrayfun(@num2str, scen_nx,  'UniformOutput',false), ',') ']'];
ny_js       = ['[' strjoin(arrayfun(@num2str, scen_ny,  'UniformOutput',false), ',') ']'];
dist_js     = ['[' strjoin(arrayfun(@(v)sprintf('%.3f',v), scen_dkm, 'UniformOutput',false), ',') ']'];
stat_gt_js  = ['[' strjoin(arrayfun(@(v)sprintf('%.4f',v), static_Gt, 'UniformOutput',false), ',') ']'];
stat_lbl_js = ['["' strjoin(static_lbls, '","') '"]'];
stat_col_js = ['["' strjoin(static_cols, '","') '"]'];
stat_dc_js  = ['[' strjoin(arrayfun(@(v)sprintf('%.3f',v), static_dc, 'UniformOutput',false), ',') ']'];

if any_uncertain
    mode_str = 'Det + Joint MC (P10/P50/P90)';
else
    mode_str = 'Deterministic';
end

chart_title = sprintf('Dynamic vs Static Storage Capacity: %s', mode_str);
chart_sub   = sprintf('%s  |  best pattern & spacing per well count  |  %d well counts  |  %d static ceilings', ...
    char(OUT.center_point_type), N, n_static);

%% 6. Write HTML
fid = fopen(html_out, 'w', 'n', 'UTF-8');
if fid < 0, error('Cannot open "%s" for writing.', html_out); end
w = @(s) fprintf(fid, '%s\n', s);
fprintf(fid, '<!DOCTYPE html><html lang="en"><head>\n');
fprintf(fid, '<meta charset="UTF-8">\n');
fprintf(fid, '<meta name="viewport" content="width=device-width,initial-scale=1.0">\n');
fprintf(fid, '<title>%s</title>\n', chart_title);
fprintf(fid, '<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>\n');
w('<style>');
w('*{box-sizing:border-box;margin:0;padding:0}');
w('body{font-family:"Segoe UI",system-ui,sans-serif;background:#eef1f7;min-height:100vh;');
w('  display:flex;flex-direction:column;align-items:center;padding:24px 16px 40px}');
w('.card{background:#fff;border-radius:14px;box-shadow:0 2px 20px rgba(0,0,0,.10);');
w('  padding:28px 32px 20px;width:100%;max-width:1380px}');
w('h1{font-size:1.25rem;font-weight:700;color:#1a2233;letter-spacing:-.01em}');
w('.sub{font-size:.80rem;color:#6b7a99;margin-top:4px;margin-bottom:18px}');
w('.wrap{position:relative;width:100%;height:540px}');
w('.hint{text-align:center;font-size:.74rem;color:#8899bb;margin-top:10px}');
w('#tt{position:fixed;pointer-events:none;background:rgba(14,22,42,.96);color:#fff;');
w('  border-radius:10px;padding:13px 16px;font-size:.79rem;line-height:1.8;');
w('  min-width:260px;box-shadow:0 8px 28px rgba(0,0,0,.32);z-index:9999;display:none;');
w('  backdrop-filter:blur(6px);border:1px solid rgba(255,255,255,.10)}');
w('.tt-hd{font-size:.88rem;font-weight:700;color:#e2eaf8;');
w('  border-bottom:1px solid rgba(255,255,255,.14);padding-bottom:6px;margin-bottom:6px}');
w('.tt-meta{color:#8aacce;font-size:.74rem;margin-bottom:8px}');
w('.tt-sec{font-size:.70rem;color:#6a9ab8;text-transform:uppercase;letter-spacing:.06em;');
w('  margin:8px 0 3px;padding-top:4px;border-top:1px solid rgba(255,255,255,.07)}');
w('.tt-row{display:flex;align-items:center;gap:8px;margin:2px 0}');
w('.tt-dot{width:9px;height:9px;border-radius:50%;flex-shrink:0;border:1.5px solid rgba(255,255,255,.5)}');
w('.tt-nm{flex:1;color:#b8cfe8;font-size:.77rem}');
w('.tt-vl{font-weight:700;color:#fff;font-size:.81rem}');
w('.tt-ratio{font-size:.74rem;color:#8fc8e0;margin-left:3px}');
fprintf(fid, '</style></head><body><div class="card">\n');
fprintf(fid, '<h1>%s</h1>\n', chart_title);
fprintf(fid, '<div class="sub">%s</div>\n', chart_sub);
fprintf(fid, '<div class="wrap"><canvas id="ch"></canvas></div>\n');
w('<p class="hint">&#9432; Hover over a well count:');
w('  solid/dashed = dynamic capacity (varies per well count) &bull;');
w('  dotted = static ceiling (constant, from formation DeltaPcr) &bull;');
w('  ratio = dynamic / static x 100%</p>');
fprintf(fid, '</div><div id="tt"></div><script>\n');
fprintf(fid, 'const mW   = %s;\n', wells_js);
fprintf(fid, 'const mNX  = %s;\n', nx_js);
fprintf(fid, 'const mNY  = %s;\n', ny_js);
fprintf(fid, 'const mD   = %s;\n', dist_js);
fprintf(fid, 'const mSGt = %s;\n', stat_gt_js);
fprintf(fid, 'const mSL  = %s;\n', stat_lbl_js);
fprintf(fid, 'const mSC  = %s;\n', stat_col_js);
fprintf(fid, 'const mSDc = %s;\n', stat_dc_js);
fprintf(fid, 'const nDyn          = %d;\n', n_dyn);
fprintf(fid, 'const nStatic       = %d;\n', n_static);
fprintf(fid, 'const N             = %d;\n', N);
fprintf(fid, 'const permOnlyStatic= %s;\n', ifelse(perm_only_static,'true','false'));
fprintf(fid, 'const ctx = document.getElementById("ch").getContext("2d");\n');
fprintf(fid, 'const chart = new Chart(ctx, {\n');
fprintf(fid, '  type: "line",\n');
fprintf(fid, '  data: { labels: [%s], datasets: [%s] },\n', xlabs, datasets_js);
fprintf(fid, '  options: {\n');
fprintf(fid, '    responsive: true, maintainAspectRatio: false,\n');
fprintf(fid, '    interaction: { mode: "index", intersect: false },\n');
fprintf(fid, '    plugins: {\n');
fprintf(fid, '      legend: { display: true, position: "top", align: "end",\n');
fprintf(fid, '        labels: { boxWidth: 28, boxHeight: 3, padding: 11,\n');
fprintf(fid, '          font: { size: 10, weight: "500" }, color: "#2c3a55",\n');
fprintf(fid, '          usePointStyle: true, pointStyle: "circle" } },\n');
fprintf(fid, '      tooltip: { enabled: false } },\n');
fprintf(fid, '    scales: {\n');
fprintf(fid, '      x: { title: { display: true, text: "Number of Wells",\n');
fprintf(fid, '             font: { size: 12, weight: "600" }, color: "#334" },\n');
fprintf(fid, '           ticks: { color: "#556", font: { size: 11 },\n');
fprintf(fid, '             maxRotation: %d, autoSkip: true, maxTicksLimit: %d },\n', ...
    ifelse_num(N>30,55,0), min(N,60));
fprintf(fid, '           grid:  { color: "rgba(100,120,160,.11)" } },\n');
fprintf(fid, '      y: { type: "logarithmic",\n');     % log y-axis
fprintf(fid, '           title: { display: true, text: "CO2 Storage Capacity [Gt]",\n');
fprintf(fid, '             font: { size: 12, weight: "600" }, color: "#334" },\n');
fprintf(fid, '           ticks: { color: "#556", font: { size: 11 } },\n');
fprintf(fid, '           grid:  { color: "rgba(100,120,160,.11)" } } },\n');
fprintf(fid, '    animation: { duration: 700, easing: "easeOutQuart" } } });\n');

% custom tooltip
fprintf(fid, 'const tt = document.getElementById("tt");\n');
fprintf(fid, 'document.getElementById("ch").addEventListener("mousemove", function(e) {\n');
fprintf(fid, '  const pts = chart.getElementsAtEventForMode(e,"index",{intersect:false},true);\n');
fprintf(fid, '  if (!pts.length) { tt.style.display="none"; return; }\n');
fprintf(fid, '  const i = pts[0].index;\n');
fprintf(fid, '  const ww = mW[i], d = mD[i].toFixed(2), pat = mNX[i]+"x"+mNY[i];\n');
fprintf(fid, '  let dynRows = "", statRows = "";\n');
fprintf(fid, '  // dynamic lines\n');
fprintf(fid, '  for (let li = 0; li < nDyn; li++) {\n');
fprintf(fid, '    const ds = chart.data.datasets[li], v = ds.data[i];\n');
fprintf(fid, '    if (v === null || v === undefined) continue;\n');
fprintf(fid, '    dynRows += `<div class="tt-row">`\n');
fprintf(fid, '      + `<div class="tt-dot" style="background:${ds.borderColor}"></div>`\n');
fprintf(fid, '      + `<span class="tt-nm">${ds.label}</span>`\n');
fprintf(fid, '      + `<span class="tt-vl">${v.toFixed(4)} Gt</span></div>`;\n');
fprintf(fid, '  }\n');
fprintf(fid, '  // static ceilings + ratio\n');
fprintf(fid, '  if (permOnlyStatic) {\n');
fprintf(fid, '    // single shared ceiling\n');
fprintf(fid, '    const sv = mSGt[0], sl = mSL[0], sc = mSC[0], dc = mSDc[0];\n');
fprintf(fid, '    statRows += `<div class="tt-row">`\n');
fprintf(fid, '      + `<div class="tt-dot" style="background:${sc};border-radius:2px"></div>`\n');
fprintf(fid, '      + `<span class="tt-nm">${sl} - DeltaPcr=${dc} MPa</span>`\n');
fprintf(fid, '      + `<span class="tt-vl">${sv.toFixed(4)} Gt</span></div>`;\n');
fprintf(fid, '    statRows += `<div style="font-size:.72rem;color:#8aacce;margin:4px 0 2px">Utilisation (dyn / static):</div>`;\n');
fprintf(fid, '    for (let li = 0; li < nDyn; li++) {\n');
fprintf(fid, '      const ds = chart.data.datasets[li], v = ds.data[i];\n');
fprintf(fid, '      if (v === null || v === undefined) continue;\n');
fprintf(fid, '      const ratio = (v/sv*100).toFixed(1);\n');
fprintf(fid, '      statRows += `<div class="tt-row">`\n');
fprintf(fid, '        + `<div class="tt-dot" style="background:${ds.borderColor}"></div>`\n');
fprintf(fid, '        + `<span class="tt-nm">${ds.label}</span>`\n');
fprintf(fid, '        + `<span class="tt-ratio">${ratio}%%</span></div>`;\n');
fprintf(fid, '    }\n');
fprintf(fid, '  } else {\n');
fprintf(fid, '    // multiple ceilings: match to dynamic line by color\n');
fprintf(fid, '    for (let si = 0; si < nStatic; si++) {\n');
fprintf(fid, '      const sv = mSGt[si], sl = mSL[si], sc = mSC[si], dc = mSDc[si];\n');
fprintf(fid, '      let dynV = null;\n');
fprintf(fid, '      for (let li = 0; li < nDyn; li++) {\n');
fprintf(fid, '        if (chart.data.datasets[li].borderColor === sc) {\n');
fprintf(fid, '          dynV = chart.data.datasets[li].data[i]; break; }\n');
fprintf(fid, '      }\n');
fprintf(fid, '      const ratio = (dynV !== null && dynV !== undefined)\n');
fprintf(fid, '        ? ` (${(dynV/sv*100).toFixed(1)}%%)` : "";\n');
fprintf(fid, '      statRows += `<div class="tt-row">`\n');
fprintf(fid, '        + `<div class="tt-dot" style="background:${sc}"></div>`\n');
fprintf(fid, '        + `<span class="tt-nm">${sl} DeltaPcr=${dc} MPa</span>`\n');
fprintf(fid, '        + `<span class="tt-vl">${sv.toFixed(4)} Gt</span>`\n');
fprintf(fid, '        + `<span class="tt-ratio">${ratio}</span></div>`;\n');
fprintf(fid, '    }\n');
fprintf(fid, '  }\n');
fprintf(fid, '  tt.innerHTML = `<div class="tt-hd">${ww} wells (best scenario)</div>`\n');
fprintf(fid, '    + `<div class="tt-meta">Pattern: ${pat} - Spacing: ${d} km</div>`\n');
fprintf(fid, '    + `<div class="tt-sec">Dynamic capacity</div>${dynRows}`\n');
fprintf(fid, '    + `<div class="tt-sec">Static ceiling (dyn/static ratio)</div>${statRows}`;\n');
fprintf(fid, '  tt.style.display = "block";\n');
fprintf(fid, '  const gap=18, tw=tt.offsetWidth, th=tt.offsetHeight;\n');
fprintf(fid, '  let tx=e.clientX+gap, ty=e.clientY+gap;\n');
fprintf(fid, '  if (tx+tw > window.innerWidth-8)  tx = e.clientX-tw-gap;\n');
fprintf(fid, '  if (ty+th > window.innerHeight-8) ty = e.clientY-th-gap;\n');
fprintf(fid, '  tt.style.left = tx+"px"; tt.style.top = ty+"px";\n');
fprintf(fid, '});\n');
fprintf(fid, 'document.getElementById("ch").addEventListener("mouseleave",()=>{\n');
fprintf(fid, '  tt.style.display="none";\n');
fprintf(fid, '});\n');
fprintf(fid, '</script></body></html>\n');
fclose(fid);
fprintf('Saved: %s\n', html_out);
fprintf('  %d well counts | %d dynamic lines | %d static ceiling lines\n', N, n_dyn, n_static);
try, web(html_out, '-browser'); catch, fprintf('Open manually: %s\n', html_out); end
end

%% Helpers
function v = safe_field(s, fname, default)
    if isfield(s, fname), v = s.(fname); else, v = default; end
end

function v = safe_scalar(s, fname, default)
    if isfield(s, fname) && isscalar(s.(fname)) && isfinite(s.(fname)) && s.(fname) > 0
        v = s.(fname);
    else
        v = default;
    end
end

function out = ifelse(cond, a, b)
    if cond, out = a; else, out = b; end
end

function out = ifelse_num(cond, a, b)
    if cond, out = a; else, out = b; end
end