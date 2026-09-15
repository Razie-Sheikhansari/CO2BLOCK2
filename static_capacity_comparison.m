% static_capacity_comparison.m
% Static (compressibility-based) CO2 storage capacity vs CO2BLOCK dynamic results.
% Requires CO2BLOCK_results.mat. Run CO2BLOCK.m first.

clearvars; close all;

%% load
if ~isfile('CO2BLOCK_results.mat')
    error('CO2BLOCK_results.mat not found. Run CO2BLOCK.m first.');
end
load('CO2BLOCK_results.mat', 'OUT', 'S', 'R');

if R.rc == inf
    error('domain_type is OPEN. Static formula requires a closed boundary.');
end

%% determine uncertainty mode
is_fault = (OUT.center_point_type == "Fault");

stress_active = isfield(OUT,'stress_uncertain')  && OUT.stress_uncertain  == "yes";
fault_active  = isfield(OUT,'fault_uncertainty') && OUT.fault_uncertainty == "yes";
perm_active   = isfield(OUT,'perm_uncertain')    && logical(OUT.perm_uncertain);
is_det        = ~stress_active && ~fault_active && ~perm_active;

if is_det
    n_modes    = 1;
    mode_label = {'Deterministic'};
    mode_short = {'Det'};
else
    n_modes    = 3;
    mode_label = {'P90  (conservative)', 'P50  (median)', 'P10  (optimistic)'};
    mode_short = {'P90', 'P50', 'P10'};
end

%% reservoir parameters
area_res_m2 = R.area_res * 1e6;   % km2 -> m2
thick       = R.thick;
por         = R.por;
compr       = R.compr;             % Pa-1 (includes porosity: compr = (cr+cw)*por)
dens_c      = R.dens_c;           % kg/m3
time_yr     = S.time_yr;
use_sf      = strcmp(R.use_shape_factor, 'yes');

V_bulk_m3 = area_res_m2 * thick;

rc_eff_m = R.rc;
if use_sf, rc_eff_m = sqrt(R.x_length * R.y_length / pi); end

%% Dp_cr per mode [MPa]
% Deterministic: use OUT.deltap_cr_det (minimum across faults in Fault mode).
% MC modes: compute percentiles from OUT.deltap_cr_samples;
%   Fault mode samples are a cell {nr_faults x 1}, each [n_MC x 1] -- take min per realization.
deltap_cr_MPa = NaN(n_modes, 1);

if is_det
    dc = OUT.deltap_cr_det;
    dc = dc(dc > 0);
    if ~isempty(dc), deltap_cr_MPa(1) = min(dc, [], 'omitnan'); end
else
    samps = OUT.deltap_cr_samples;
    if iscell(samps)
        % fault mode: each cell is one fault's [n_MC x 1] samples
        samp_mat    = cell2mat(samps(:)');
        dc_per_real = min(samp_mat, [], 2, 'omitnan');   % min across faults per realization
    else
        dc_per_real = samps(:);
    end
    dc_per_real(dc_per_real <= 0) = NaN;

    pcts = prctile(dc_per_real, [90, 50, 10]);
    deltap_cr_MPa(1) = pcts(1);   % P90
    deltap_cr_MPa(2) = pcts(2);   % P50
    deltap_cr_MPa(3) = pcts(3);   % P10
end

deltap_cr_Pa = deltap_cr_MPa * 1e6;

%% static capacity [Gt]
% M_static = V_bulk * compr * Dp_cr * dens_c  (por already embedded in compr)
M_static_Gt = V_bulk_m3 * compr * deltap_cr_Pa * dens_c / 1e12;

%% dynamic capacity [Gt]
well_list   = OUT.well_list;
d_list      = OUT.d_list;
x_grid_list = OUT.x_grid_list;
y_grid_list = OUT.y_grid_list;

% field names per mode
if is_det
    res_fields   = {'V_M_det'};
    fault_fields = {'V_M_fault_det'};
else
    res_fields   = {'V_M_joint_P90', 'V_M_joint_P50', 'V_M_joint_P10'};
    fault_fields = {'V_M_joint_fault_P90', 'V_M_joint_fault_P50', 'V_M_joint_fault_P10'};
end

M_dyn_res_Gt   = NaN(n_modes, 1);
M_dyn_fault_Gt = NaN(n_modes, 1);
best_nwells    = NaN(n_modes, 1);
best_dist_km   = NaN(n_modes, 1);
best_nx        = NaN(n_modes, 1);
best_ny        = NaN(n_modes, 1);

for k = 1:n_modes
    if isfield(OUT, res_fields{k})
        V = OUT.(res_fields{k});
        if ~isempty(V) && any(isfinite(V(:)))
            [M_dyn_res_Gt(k), lin_idx] = max(V(:), [], 'omitnan');
            [br, bc]       = ind2sub(size(V), lin_idx);
            best_nwells(k) = well_list(br);
            best_dist_km(k)= d_list(bc);
            best_nx(k)     = x_grid_list(br);
            best_ny(k)     = y_grid_list(br);
        end
    end
    if is_fault && isfield(OUT, fault_fields{k})
        Vf = OUT.(fault_fields{k});
        if ~isempty(Vf), M_dyn_fault_Gt(k) = max(Vf(:), [], 'omitnan'); end
    end
end

%% print summary
if is_det, unc_label = 'Deterministic'; else, unc_label = 'MC (P10/P50/P90)'; end

fprintf('\n');
fprintf('======================================================================\n');
fprintf('  STATIC vs DYNAMIC CO2 STORAGE CAPACITY COMPARISON\n');
fprintf('  %s  |  Fault mode = %s\n', unc_label, yn(is_fault));
fprintf('======================================================================\n\n');

fprintf('  area_res = %.2f km2  |  thick = %.1f m  |  por = %.4f\n', R.area_res, thick, por);
fprintf('  compr = %.3e Pa-1  |  dens_c = %.1f kg/m3  |  t = %.0f yr\n', compr, dens_c, time_yr);
fprintf('  V_bulk = %.4e m3  |  shape_factor = %s\n\n', V_bulk_m3, upper(R.use_shape_factor));

fprintf('  %-20s  %10s  %12s  %12s  %12s  %12s\n', ...
        'Mode','Dp_cr[MPa]','Static[Gt]','Dyn-Res[Gt]','Dyn/Stat[%]','Dyn-Fault[Gt]');
fprintf('  %s\n', repmat('-',1,85));

for k = 1:n_modes
    ratio_pct = M_dyn_res_Gt(k) / M_static_Gt(k) * 100;
    fault_str = sprintf('%12s', 'n/a');
    if is_fault && ~isnan(M_dyn_fault_Gt(k))
        fault_str = sprintf('%12.4f', M_dyn_fault_Gt(k));
    end
    flag = ''; if ratio_pct > 100, flag = '  <- (*)'; end
    fprintf('  %-20s  %10.3f  %12.4f  %12.4f  %12.1f  %s%s\n', ...
        mode_label{k}, deltap_cr_MPa(k), M_static_Gt(k), M_dyn_res_Gt(k), ratio_pct, fault_str, flag);
end
fprintf('\n');

time_s        = time_yr * 86400 * 365;
R_infl_km     = sqrt(2.246 * R.perm * time_s / (R.visc_w * R.compr)) / 1000;
rc_km         = rc_eff_m / 1000;
boundary_felt = R_infl_km >= rc_km;

if any(M_dyn_res_Gt > M_static_Gt)
    fprintf('  WARNING (*) Dynamic > Static\n');
    if ~boundary_felt
        fprintf('     Boundary not yet felt (R_infl=%.1f km < rc=%.1f km).\n', R_infl_km, rc_km);
        fprintf('     Reservoir behaves as open in %.0f yr. Static is the true long-term limit.\n\n', time_yr);
    else
        fprintf('     Boundary IS felt (R_infl=%.1f km >= rc=%.1f km).\n', R_infl_km, rc_km);
        fprintf('     Check maxQ, perm, and compr values.\n\n');
    end
end

fprintf('  Best scenario (reservoir-limited):\n');
for k = 1:n_modes
    fprintf('    %s: %dx%d grid (%d wells) | %.1f km spacing\n', ...
        mode_short{k}, best_nx(k), best_ny(k), best_nwells(k), best_dist_km(k));
end

fprintf('\n  R_influence = %.1f km  |  rc = %.1f km  |  Boundary felt = %s\n', R_infl_km, rc_km, yn(boundary_felt));
max_ratio = max(M_dyn_res_Gt ./ M_static_Gt, [], 'omitnan') * 100;
if ~boundary_felt
    fprintf('  REGIME: BOUNDARY NOT YET FELT -- extend time_yr to converge to static limit.\n');
elseif max_ratio < 80
    fprintf('  REGIME: PERMEABILITY-LIMITED (%.0f%%) -- more wells or time would help.\n', max_ratio);
elseif max_ratio <= 110
    fprintf('  REGIME: STORAGE-LIMITED (%.0f%%) -- near compressibility bound.\n', max_ratio);
else
    fprintf('  REGIME: OVER-INJECTION WARNING (%.0f%%) -- review maxQ/perm/compr.\n', max_ratio);
end
fprintf('======================================================================\n\n');

%% bar chart
figure('Name','Static vs Dynamic Storage Capacity','NumberTitle','off','Color','w');

if is_fault && any(~isnan(M_dyn_fault_Gt))
    bar_data    = [M_static_Gt, M_dyn_res_Gt, M_dyn_fault_Gt];
    leg_entries = {'Static bound','Dynamic - reservoir','Dynamic - fault-constrained'};
    bar_colors  = [0.75 0.88 1.0; 0.20 0.55 0.85; 0.85 0.35 0.20];
else
    bar_data    = [M_static_Gt, M_dyn_res_Gt];
    leg_entries = {'Static bound','Dynamic - reservoir'};
    bar_colors  = [0.75 0.88 1.0; 0.20 0.55 0.85];
end

n_bars        = size(bar_data, 2);
bar_data_plot = bar_data;
if n_modes == 1, bar_data_plot = [bar_data; NaN(1, n_bars)]; end

b_h = bar(bar_data_plot, 0.65);
for bi = 1:n_bars
    b_h(bi).FaceColor = bar_colors(bi,:);
    b_h(bi).EdgeColor = 'none';
end

hold on;
x_centers = get(b_h(1), 'XData');
offsets   = linspace(-(n_bars-1)/2, (n_bars-1)/2, n_bars) .* (0.65/n_bars);
for gi = 1:n_modes
    for bi = 1:n_bars
        val = bar_data(gi, bi);
        if isnan(val), continue; end
        text(x_centers(gi) + offsets(bi), val + 0.01*max(M_static_Gt,[],'omitnan'), ...
             sprintf('%.3f',val), 'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
    end
end

xticks(1:n_modes); xticklabels(mode_short);
xlabel(sprintf('Mode  (%s)', unc_label));
ylabel('CO2 Storage Capacity  [Gt]');
title(sprintf('Static vs Dynamic Storage Capacity  |  t = %.0f yr', time_yr));
legend(leg_entries,'Location','northeast');
grid on; box on;
y_top = max(M_static_Gt,[],'omitnan') * 1.30;
if isnan(y_top) || y_top == 0, y_top = 1; end
ylim([0, y_top]); set(gca,'FontSize',11);
hold off;

%% helper
function s = yn(flag)
    if (islogical(flag)||isnumeric(flag)), s = 'NO'; if flag, s = 'YES'; end
    elseif strcmpi(char(flag),'yes'),      s = 'YES';
    else,                                  s = 'NO';
    end
end