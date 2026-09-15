% GEOS_QMAX_COMPARISON
% Compares Q_max and V_max from a GEOS pressure map against the
% CO2BLOCK2 deterministic result for the same well/distance scenario.
% Run CO2BLOCK2.m first, then run this script (F5 or Run button).

% -------------------------------------------------------------------------
% USER INPUTS: edit these to match your GEOS run
n_wells          = 4;        % number of injection wells in GEOS
d_geos_km        = 0.514;    % inter-well distance [km]
delta_p_max_geos = 5.67;     % max overpressure from GEOS map [MPa]
% -------------------------------------------------------------------------

%% load results if not already in workspace
if ~exist('OUT','var') || ~exist('S','var') || ~exist('R','var')
    if exist('CO2BLOCK2_results.mat','file')
        fprintf('Loading CO2BLOCK_results.mat ...\n');
        load('CO2BLOCK2_results.mat', 'OUT', 'S', 'R');
    else
        error(['OUT, S, R not found in workspace and CO2BLOCK_results.mat not found.\n' ...
               'Run CO2BLOCK2.m first, then run this script.']);
    end
end

%% read settings
time_yr   = S.time_yr;
well_list = OUT.well_list(:);
d_list    = OUT.d_list(:);   % [km]

% critical overpressure (deterministic)
if isfield(OUT,'deltap_cr_ref') && isfinite(OUT.deltap_cr_ref) && OUT.deltap_cr_ref > 0
    deltap_cr = OUT.deltap_cr_ref;   % [MPa]
else
    error(['OUT.deltap_cr_ref not found or invalid.\n' ...
           'Make sure CO2BLOCK2 was run with a valid geomechanical setup.']);
end

%% find matching well scenario
wid = find(well_list == n_wells, 1);
if isempty(wid)
    fprintf('\nERROR: n_wells = %d not found in OUT.well_list.\n', n_wells);
    fprintf('  Available well counts: %s\n', num2str(well_list'));
    return
end

%% find matching distance scenario
[dist_err, did] = min(abs(d_list - d_geos_km));
d_matched_km    = d_list(did);
if dist_err > 0.05
    fprintf('WARNING: exact distance %.3f km not found. Using nearest: %.3f km\n', ...
        d_geos_km, d_matched_km);
end

%% reference pressure and scaling
% p_sup_vec is the superposed overpressure at the reference point, at Q0
p_sup_co2block = OUT.p_sup_vec(wid, did);   % [MPa]

Q0_m3s  = OUT.Q0_vec(wid);                          % [m3/s per well]
Q0_Mtyr = Q0_m3s * 86400 * 365 * R.dens_c / 1e9;   % [Mt/yr per well]

scale_factor_geos   = deltap_cr / delta_p_max_geos;
scale_factor_co2blk = deltap_cr / p_sup_co2block;   %#ok<NASGU>

Q_max_geos_tyr = Q0_Mtyr * 1e6 * scale_factor_geos;
V_max_geos_t   = Q_max_geos_tyr * n_wells * time_yr;

Q_co2block_tyr = OUT.Q_M_each_1(wid, did) * 1e6;
V_co2block_t   = OUT.V_M_1(wid, did)      * 1e9;

%% print domain properties
sep = repmat('-', 1, 62);

fprintf('\n%s\n', sep);
fprintf('  Domain & Reservoir Properties\n');
fprintf('%s\n', sep);
if isfield(R,'x_length') && isfield(R,'y_length')
    fprintf('  Domain size    : %.0f m  x  %.0f m\n', R.x_length, R.y_length);
end
fprintf('  Domain area    : %.4f km2\n',   R.area_res);
fprintf('  Thickness      : %.1f m\n',     R.thick);
fprintf('  Porosity       : %.4f\n',       R.por);
fprintf('  Permeability   : %.4f mD\n',    R.perm / 9.869e-16);
fprintf('  CO2 density    : %.2f kg/m3\n', R.dens_c);
fprintf('  Injection time : %.4f yr\n',    time_yr);
fprintf('%s\n', sep);

%% print comparison table
fprintf('  GEOS vs CO2BLOCK2 Comparison\n');
fprintf('  Wells : %d  |  Distance : %.3f km  |  Time : %.4f yr\n', ...
    n_wells, d_matched_km, time_yr);
fprintf('%s\n', sep);
fprintf('  %-22s  %20s  %16s\n', 'Source', 'Q_max/well [t/yr]', 'V_max [t]');
fprintf('  %-22s  %20s  %16s\n', repmat('-',1,22), repmat('-',1,20), repmat('-',1,16));
fprintf('  %-22s  %20.2f  %16.2f\n', 'GEOS-derived',   Q_max_geos_tyr, V_max_geos_t);
fprintf('  %-22s  %20.2f  %16.2f\n', 'CO2BLOCK2 (det.)', Q_co2block_tyr, V_co2block_t);
fprintf('%s\n', sep);

if Q_max_geos_tyr > 0 && isfinite(Q_max_geos_tyr) && isfinite(Q_co2block_tyr)
    dQ_pct = (Q_co2block_tyr - Q_max_geos_tyr) / Q_max_geos_tyr * 100;
    dV_pct = (V_co2block_t   - V_max_geos_t)   / V_max_geos_t   * 100;
    fprintf('  %-22s  %+19.4f%%  %+15.4f%%\n', 'CO2BLOCK2 vs GEOS', dQ_pct, dV_pct);
    fprintf('%s\n', sep);
end

%% diagnostics
fprintf('  Diagnostics:\n');
fprintf('    Q0 reference (%d wells)       : %.6f Mt/yr/well\n', n_wells, Q0_Mtyr);
fprintf('    Q0 reference (%d wells)       : %.4f t/yr/well\n',  n_wells, Q0_Mtyr*1e6);
fprintf('    delta_P_cr   (CO2BLOCK2)       : %.4f MPa\n', deltap_cr);
fprintf('    delta_P_max  (GEOS)           : %.4f MPa\n', delta_p_max_geos);
fprintf('    delta_P_sup  (CO2BLOCK2 at Q0) : %.4f MPa\n', p_sup_co2block);

if ~isfinite(Q_co2block_tyr) || isnan(Q_co2block_tyr)
    fprintf('\n  NOTE: CO2BLOCK2 result is NaN for this scenario.\n');
    fprintf('  Distance %.3f km may exceed d_max for %d wells (%.3f km).\n', ...
        d_matched_km, n_wells, OUT.d_max(wid));
end

fprintf('\n');