function G = read_geomech()
% Reads geomechanical parameters from Geomech_parameters.xlsx
% and fault geometry from Faults.xlsx.
%
% Each stress component (pres, Sv, Shmin, SHmax) can be given as:
%   - a gradient (*_grad_seismic) plus optional offset, or
%   - a single known value (*_known), with type 'effective' or 'total'
% If both are provided, the direct value wins.
%
% depth_mean, depth_seismic, and SHmax_dir_seismic are always required.
% Shmin and SHmax must always be specified (no fallback default).
%
% T_mean_known and pres_mean_known are optional; if given, they bypass
% the gradient calculation in compute_fluid_properties.
% depth_water defaults to 0 (onshore) if not set.
%
% Note: in direct mode, a_*_grad is interpreted as +/- MPa (not MPa/km)
% and is back-converted here so the Monte Carlo code stays unchanged.

%% read file into parameter map
raw    = readcell('Geomech_parameters.xlsx', 'Sheet', 1, 'Range', 'A:B');
raw    = raw(~cellfun(@(x) ~ischar(x) || isempty(strtrim(x)), raw(:,1)), :);
params = containers.Map(strtrim(raw(:,1)), raw(:,2));

%% helper functions
function v = num(key),  v = double(params(key));          end
function v = str(key),  v = strtrim(char(params(key)));   end
function v = num_opt(key, def)
    if has_value(key), v = double(params(key)); else, v = def; end
end
function v = str_opt(key, def)
    if has_value(key), v = strtrim(char(params(key))); else, v = def; end
end
function tf = has_value(key)
    tf = isKey(params, key) && ~cell_is_blank(params(key));
end

%% Reservoir geometry
G.depth_mean  = num('depth_mean');         % required, no default
G.depth       = num_opt('depth', 0);       % optional; handled in compute_fluid_properties if 0
G.depth_water = num_opt('depth_water', 0); % 0 for onshore

if ~has_value('depth_water')
    fprintf('read_geomech: depth_water not set, defaulting to 0 m (onshore).\n');
    fprintf('  For offshore: fill depth_water, or supply T_mean_known + pres_mean_known.\n');
end

%% Fluid property shortcuts (optional)
% If set, compute_fluid_properties uses these directly and skips the gradient path.
G.T_mean_known    = num_opt('T_mean_known',    NaN);  % [C]
G.pres_mean_known = num_opt('pres_mean_known', NaN);  % [MPa]

if isfinite(G.T_mean_known)
    fprintf('read_geomech: T_mean_known = %.1f C (skipping temp gradient for fluid props)\n', ...
        G.T_mean_known);
end
if isfinite(G.pres_mean_known)
    fprintf('read_geomech: pres_mean_known = %.2f MPa (skipping pres gradient for fluid props)\n', ...
        G.pres_mean_known);
end

%% Pressure and thermal gradients
% Used by compute_fluid_properties when T_mean_known / pres_mean_known are not set.
G.pres_grad   = num_opt('pres_grad',  10.2) / 1e3;  % MPa/km -> MPa/m
G.temp_surf   = num_opt('temp_surf',  15);           % [C]
G.temp_grad   = num_opt('temp_grad',  33);           % [C/km]
G.pres_offset = num_opt('pres_offset', 0);           % [MPa] gradient intercept

%% Always-required stress fields
G.depth_seismic     = num('depth_seismic');
G.SHmax_dir_seismic = num_opt('SHmax_dir_seismic', NaN);

depth_eff = G.depth_seismic - G.depth_water;

%% Resolve pore pressure (must come first; needed for effective-to-total conversion)
direct_p = has_value('pres_known');
grad_p   = has_value('pres_grad_seismic');

if direct_p
    p0_val = double(params('pres_known'));
    G.pres_offset       = p0_val;
    G.pres_grad_seismic = 0;
    if grad_p
        fprintf('read_geomech: WARNING: both pres_known and pres_grad_seismic set; using direct value\n');
    end
    fprintf('read_geomech: pressure: direct value %.2f MPa at depth %.0f m\n', ...
        p0_val, G.depth_seismic);
elseif grad_p
    G.pres_grad_seismic = num('pres_grad_seismic') / 1e3;
else
    % no pressure input: fall back to hydrostatic
    G.pres_grad_seismic = 10.2 / 1e3;
    fprintf('read_geomech: WARNING: no pressure data, defaulting to hydrostatic 10.2 MPa/km\n');
end

p0 = G.pres_offset + G.pres_grad_seismic * depth_eff;
fprintf('read_geomech: pore pressure at reference depth = %.2f MPa\n', p0);

%% Resolve Sv
direct_Sv = has_value('Sv_known');
grad_Sv   = has_value('Sv_grad_seismic');

if direct_Sv
    Sv_val  = double(params('Sv_known'));
    Sv_type = str_opt('Sv_known_type', 'total');
    if ~has_value('Sv_known_type')
        fprintf('read_geomech: WARNING: Sv_known_type not set, assuming "total"\n');
    end
    if strcmpi(Sv_type, 'effective')
        Sv_val = Sv_val + p0;
        fprintf('read_geomech: Sv: effective %.2f MPa + p0 %.2f = total %.2f MPa\n', ...
            Sv_val - p0, p0, Sv_val);
    else
        fprintf('read_geomech: Sv: total value %.2f MPa\n', Sv_val);
    end
    G.Sv_offset       = Sv_val;
    G.Sv_grad_seismic = 0;
    if grad_Sv
        fprintf('read_geomech: WARNING: both Sv_known and Sv_grad_seismic set; using direct value\n');
    end
elseif grad_Sv
    G.Sv_offset       = num_opt('Sv_offset', 0);
    G.Sv_grad_seismic = num('Sv_grad_seismic') / 1e3;
else
    % no Sv input: fall back to typical overburden gradient
    G.Sv_offset       = num_opt('Sv_offset', 0);
    G.Sv_grad_seismic = 23 / 1e3;
    fprintf('read_geomech: WARNING: no Sv data, defaulting to overburden 23 MPa/km\n');
end

%% Resolve Shmin
direct_Shmin = has_value('Shmin_known');
grad_Shmin   = has_value('Shmin_grad_seismic');

if direct_Shmin
    Sh_val  = double(params('Shmin_known'));
    Sh_type = str_opt('Shmin_known_type', 'total');
    if ~has_value('Shmin_known_type')
        fprintf('read_geomech: WARNING: Shmin_known_type not set, assuming "total"\n');
    end
    if strcmpi(Sh_type, 'effective')
        Sh_val = Sh_val + p0;
        fprintf('read_geomech: Shmin: effective %.2f MPa + p0 %.2f = total %.2f MPa\n', ...
            Sh_val - p0, p0, Sh_val);
    else
        fprintf('read_geomech: Shmin: total value %.2f MPa\n', Sh_val);
    end
    G.Shmin_offset       = Sh_val;
    G.Shmin_grad_seismic = 0;
    if grad_Shmin
        fprintf('read_geomech: WARNING: both Shmin_known and Shmin_grad_seismic set; using direct value\n');
    end
elseif grad_Shmin
    G.Shmin_offset       = num_opt('Shmin_offset', 0);
    G.Shmin_grad_seismic = num('Shmin_grad_seismic') / 1e3;
else
    error(['read_geomech: Shmin has no input. ' ...
        'Provide Shmin_grad_seismic or Shmin_known. ' ...
        'No physical default exists for horizontal stress.']);
end

%% Resolve SHmax
direct_SHmax = has_value('SHmax_known');
grad_SHmax   = has_value('SHmax_grad_seismic');

if direct_SHmax
    SH_val  = double(params('SHmax_known'));
    SH_type = str_opt('SHmax_known_type', 'total');
    if ~has_value('SHmax_known_type')
        fprintf('read_geomech: WARNING: SHmax_known_type not set, assuming "total"\n');
    end
    if strcmpi(SH_type, 'effective')
        SH_val = SH_val + p0;
        fprintf('read_geomech: SHmax: effective %.2f MPa + p0 %.2f = total %.2f MPa\n', ...
            SH_val - p0, p0, SH_val);
    else
        fprintf('read_geomech: SHmax: total value %.2f MPa\n', SH_val);
    end
    G.SHmax_offset       = SH_val;
    G.SHmax_grad_seismic = 0;
    if grad_SHmax
        fprintf('read_geomech: WARNING: both SHmax_known and SHmax_grad_seismic set; using direct value\n');
    end
elseif grad_SHmax
    G.SHmax_offset       = num_opt('SHmax_offset', 0);
    G.SHmax_grad_seismic = num('SHmax_grad_seismic') / 1e3;
else
    error(['read_geomech: SHmax has no input. ' ...
        'Provide SHmax_grad_seismic or SHmax_known. ' ...
        'No physical default exists for horizontal stress.']);
end

%% Rock strength
% rock_tens_strength is optional; defaults to rock_cohesion / 2
G.rock_cohesion      = num('rock_cohesion');
G.rock_tens_strength = num_opt('rock_tens_strength', G.rock_cohesion / 2);
G.fault_friction     = num('fault_friction');

if ~has_value('rock_tens_strength')
    fprintf('read_geomech: rock_tens_strength not set, defaulting to rock_cohesion/2 = %.2f MPa\n', ...
        G.rock_tens_strength);
end

%% Distribution types
G.f_dist_p     = str('f_dist_p');
G.f_dist_Sv    = str('f_dist_Sv');
G.f_dist_Shmin = str('f_dist_Shmin');
G.f_dist_SHmax = str('f_dist_SHmax');
G.f_dist_dir   = str('f_dist_dir');
G.f_dist_dip   = str('f_dist_dip');
G.f_dist_azi   = str('f_dist_azi');
G.f_dist_mu    = str('f_dist_mu');

%% Uncertainty magnitudes
% Gradient mode: a_*_grad is in MPa/km. Downstream MC computes:
%   half_range = depth_eff * (a_*_grad / 1000)
%
% Direct mode: user enters a_*_grad as +/- absolute MPa.
% Back-convert so the downstream formula still gives the right half-range:
%   stored_grad = user_MPa * 1000 / depth_eff
G.a_p_grad     = num('a_p_grad');
G.a_Sv_grad    = num('a_Sv_grad');
G.a_Shmin_grad = num('a_Shmin_grad');
G.a_SHmax_grad = num('a_SHmax_grad');
G.a_SHmax_dir  = num('a_SHmax_dir');
G.a_fault_dip  = num('a_fault_dip');
G.a_fault_azi  = num('a_fault_azi');
G.a_fault_mu   = num('a_fault_mu');

% Back-convert direct-mode uncertainties from absolute MPa to MPa/km equivalent
if depth_eff > 0
    if direct_p
        G.a_p_grad = G.a_p_grad * 1000 / depth_eff;
        fprintf('read_geomech: a_p_grad interpreted as +/- %.2f MPa (direct mode)\n', ...
            G.a_p_grad * depth_eff / 1000);
    end
    if direct_Sv
        G.a_Sv_grad = G.a_Sv_grad * 1000 / depth_eff;
        fprintf('read_geomech: a_Sv_grad interpreted as +/- %.2f MPa (direct mode)\n', ...
            G.a_Sv_grad * depth_eff / 1000);
    end
    if direct_Shmin
        G.a_Shmin_grad = G.a_Shmin_grad * 1000 / depth_eff;
        fprintf('read_geomech: a_Shmin_grad interpreted as +/- %.2f MPa (direct mode)\n', ...
            G.a_Shmin_grad * depth_eff / 1000);
    end
    if direct_SHmax
        G.a_SHmax_grad = G.a_SHmax_grad * 1000 / depth_eff;
        fprintf('read_geomech: a_SHmax_grad interpreted as +/- %.2f MPa (direct mode)\n', ...
            G.a_SHmax_grad * depth_eff / 1000);
    end
end

%% Reference coordinates
G.ref_X = num_opt('ref_X', 0);
G.ref_Y = num_opt('ref_Y', 0);

%% Fault geometry from Faults.xlsx
F = readtable('Faults.xlsx');
G.fault_coord_x = F.fault_coord_x - G.ref_X;
G.fault_coord_y = F.fault_coord_y - G.ref_Y;
G.fault_length  = F.fault_length;
G.fault_azi     = F.fault_azi;
G.fault_dip     = F.fault_dip;
G.nr_fault      = height(F);

% fault_coord_z is optional; needed when fault_stress_mode = 'fault_depth'
if ismember('fault_coord_z', F.Properties.VariableNames)
    G.fault_coord_z = F.fault_coord_z;
else
    G.fault_coord_z = [];
end

G.fault_stress_mode = str_opt('fault_stress_mode', 'reference');
fprintf('read_geomech: fault_stress_mode = %s\n', G.fault_stress_mode);

end

%% local helper
function tf = cell_is_blank(val)
if ismissing(val)
    tf = true;
elseif isnumeric(val) && isnan(val)
    tf = true;
elseif ischar(val) || isstring(val)
    s = strtrim(char(string(val)));
    tf = isempty(s) || strcmp(s, '-');
else
    tf = false;
end
end