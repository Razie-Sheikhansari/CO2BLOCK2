function GEO = compute_critical_overpressure(S, G, R, MC)
% Geomechanical critical overpressure per fault (or reservoir-scale).
%
%   GEO = compute_critical_overpressure(S, G, R, MC)
%
%   Computes:
%     deltap_cr_det    - deterministic value (mean stress, fixed geometry)
%     deltap_cr_samples- raw MC samples (when any uncertainty is active)
%     deltap_cr_ref    - Central reference baseline (always)
%
%   No percentile aggregation - [Joint-MC] uses raw samples directly.

center_point_type = S.center_point_type;
stress_uncertain  = S.stress_uncertain;
fault_uncertainty = S.fault_uncertainty;
n_MC              = S.n_MC;

% default outputs
GEO.mode              = char(center_point_type);
GEO.nr_fault          = 0;
GEO.fault_coord_x     = [];
GEO.fault_coord_y     = [];
GEO.fault_azi         = [];
GEO.fault_dip         = [];
GEO.fault_length      = [];
GEO.deltap_cr_det     = [];
GEO.deltap_cr_samples = [];

% reference baseline (always computed, used for ratio plot)
ref_stress        = mean_stress_state(G);
GEO.deltap_cr_ref = central_deltap_crit(ref_stress, G);
fprintf('  Dp_cr_ref = %.2f MPa\n', GEO.deltap_cr_ref);


% Central / Desired reference point
if center_point_type == "Central" || center_point_type == "Desired"

    stress = mean_stress_state(G);
    dc = central_deltap_crit(stress, G);
    GEO.deltap_cr_det = dc;

    if stress_uncertain == "yes"
        assert(nargin >= 4 && MC.stress_active, ...
            'stress_uncertain=yes but MC.stress_active=false.');

        ss   = MC.stress_samples;
        vals = nan(n_MC, 1);
        for i = 1:n_MC
            stress.p     = ss.r_p(i);
            stress.Sv    = ss.r_Sv(i);
            stress.Shmin = ss.r_Shmin(i);
            stress.SHmax = ss.r_SHmax(i);
            stress.SHdir = ss.r_dir(i);
            vals(i) = central_deltap_crit(stress, G);
        end
        vals(vals <= 0) = NaN;
        GEO.deltap_cr_samples = vals;

        fprintf('  %s  stress MC  (n=%d, Dp_cr_det=%.2f MPa)\n', center_point_type, n_MC, dc);
    else
        GEO.deltap_cr_samples = [];
        fprintf('  %s  Dp_cr = %.2f MPa\n', center_point_type, dc);
    end
    return;
end


% Fault reference point

% Stage 1: fault geometry
switch S.fault_type
    case "Deterministic"
        FAULTS = get_fault_geometry_deterministic(G);
    case "Stochastic"
        if nargin < 3
            error('compute_critical_overpressure: R required for Stochastic fault mode.');
        end
        FAULTS = get_fault_geometry_stochastic(G, R.area_res, R.x_length, R.y_length, R.thick);
    otherwise
        if S.fault_type == "-" || S.fault_type == "N/A" || strtrim(char(S.fault_type)) == ""
            error(['Fault Type is blank but Reference Point is "Fault".\n' ...
                'Set Fault Type to "Deterministic" or "Stochastic".']);
        else
            error('Unknown Fault Type "%s".', S.fault_type);
        end
end

GEO.nr_fault      = FAULTS.nr_fault;
GEO.fault_coord_x = FAULTS.fault_coord_x;
GEO.fault_coord_y = FAULTS.fault_coord_y;
GEO.fault_azi     = FAULTS.fault_azi;
GEO.fault_dip     = FAULTS.fault_dip;
if isfield(FAULTS, 'fault_length')
    GEO.fault_length = FAULTS.fault_length;
end

% Stage 2: deterministic Dp_cr (mean stress, fixed geometry)
GEO = compute_deterministic_deltap_cr(GEO, FAULTS, G);

% Stage 3: MC samples if any uncertainty is active
any_mc = (stress_uncertain == "yes") || (fault_uncertainty == "yes");

if ~any_mc
    GEO.deltap_cr_samples = {};
    fprintf('  %d faults  [%s, det]\n', FAULTS.nr_fault, S.fault_type);
else
    if stress_uncertain == "yes"
        assert(nargin >= 4 && MC.stress_active, ...
            'stress_uncertain=yes but MC.stress_active=false.');
    end
    samples = run_fault_MC(FAULTS, G, MC, S, n_MC);
    GEO = aggregate_raw_samples(GEO, samples, FAULTS);
    fprintf('  %d faults  [%s, MC n=%d]\n', FAULTS.nr_fault, S.fault_type, n_MC);
end
end


function stress = mean_stress_state(G)
depth_eff = G.depth_seismic - G.depth_water;
stress.p     = G.pres_offset  + G.pres_grad_seismic  * depth_eff;
stress.Sv    = G.Sv_offset    + G.Sv_grad_seismic    * depth_eff;
stress.Shmin = G.Shmin_offset + G.Shmin_grad_seismic * depth_eff;
stress.SHmax = G.SHmax_offset + G.SHmax_grad_seismic * depth_eff;
stress.SHdir = G.SHmax_dir_seismic;
end

function stress = fault_stress_state(G, z_fault)
% Stress state evaluated at a specific fault depth z_fault [m from sea level]
depth_eff    = z_fault - G.depth_water;
stress.p     = G.pres_offset  + G.pres_grad_seismic  * depth_eff;
stress.Sv    = G.Sv_offset    + G.Sv_grad_seismic    * depth_eff;
stress.Shmin = G.Shmin_offset + G.Shmin_grad_seismic * depth_eff;
stress.SHmax = G.SHmax_offset + G.SHmax_grad_seismic * depth_eff;
stress.SHdir = G.SHmax_dir_seismic;
end


function dc = central_deltap_crit(stress, G)
p0            = stress.p;
phi           = atan(G.fault_friction);
cohesion      = G.rock_cohesion;
tens_strength = G.rock_tens_strength;

s1 = max([stress.Sv, stress.SHmax, stress.Shmin]) - p0;
s3 = min([stress.Sv, stress.SHmax, stress.Shmin]) - p0;

if s1 <= 0
    dc = NaN;
    return;
end

theta    = (1 - sin(phi)) / (1 + sin(phi));
fprintf('  s1=%.2f, s3=%.2f, s3/s1=%.4f, theta=%.4f\n', s1, s3, s3/s1, theta);
fprintf('  cohesion=%.2f MPa, friction_angle=%.1f deg\n', cohesion, rad2deg(phi));
fprintf('  stress ratio check: s3/s1 - theta = %.4f  (negative = already critically stressed)\n', s3/s1 - theta);
fprintf('  cohesion term: %.2f MPa\n', cohesion * cos(phi) / sin(phi));
dc_shear = (s3/s1 - theta) / (1 - theta) * s1 + ...
    cohesion * cos(phi) / sin(phi);
dc_tens  = s3 + tens_strength;
dc       = min(dc_shear, dc_tens);
if dc_shear <= dc_tens
    fprintf('  --> Failure mode: SHEAR (dc_shear=%.2f, dc_tens=%.2f MPa)\n', dc_shear, dc_tens);
else
    fprintf('  --> Failure mode: TENSILE (dc_shear=%.2f, dc_tens=%.2f MPa)\n', dc_shear, dc_tens);
end
end


function GEO = compute_deterministic_deltap_cr(GEO, FAULTS, G)
use_fault_depth = strcmp(G.fault_stress_mode, 'fault_depth') && ...
    isfield(FAULTS, 'fault_coord_z') && ...
    ~isempty(FAULTS.fault_coord_z);
if ~use_fault_depth
    stress_ref = mean_stress_state(G);
end
mu  = G.fault_friction;
tens = G.rock_tens_strength;
nr  = FAULTS.nr_fault;
dc  = nan(nr, 1);

for j = 1:nr
    if use_fault_depth
        stress = fault_stress_state(G, FAULTS.fault_coord_z(j));
    else
        stress = stress_ref;
    end
    [Sn, Tau] = stress_projection( ...
        stress.SHdir, stress.SHmax, stress.Shmin, stress.Sv, ...
        FAULTS.fault_azi(j), FAULTS.fault_dip(j));
    dc_shear  = Sn - stress.p - Tau / mu;
    dc_tens   = Sn - stress.p + tens;
    dc(j)     = min(dc_shear, dc_tens);
end

dc(dc <= 0) = NaN;
GEO.deltap_cr_det = dc;
end


function samples = run_fault_MC(FAULTS, G, MC, S, n_MC)

stress_uncertain  = S.stress_uncertain;
fault_uncertainty = S.fault_uncertainty;
nr_fault          = FAULTS.nr_fault;
tens              = G.rock_tens_strength;

mean_stress = mean_stress_state(G);

if stress_uncertain == "yes"
    ss = MC.stress_samples;
end

if fault_uncertainty == "yes"
    r_dip = nan(n_MC, nr_fault);
    r_azi = nan(n_MC, nr_fault);
    for j = 1:nr_fault
        pd_dip = build_dist_angle(G.f_dist_dip, FAULTS.fault_dip(j), G.a_fault_dip);
        pd_dip = truncate(pd_dip, 0, 90);
        pd_azi = build_dist_angle(G.f_dist_azi, FAULTS.fault_azi(j), G.a_fault_azi);
        r_dip(:,j) = random(pd_dip, n_MC, 1);
        r_azi(:,j) = mod(random(pd_azi, n_MC, 1), 360);
    end
    pd_mu = build_dist_angle(G.f_dist_mu, G.fault_friction, G.a_fault_mu);
    r_mu  = random(pd_mu, n_MC, 1);
end

deltap_cr_s = cell(1, nr_fault);
deltap_cr_t = cell(1, nr_fault);
for j = 1:nr_fault
    deltap_cr_s{j} = nan(n_MC, 1);
    deltap_cr_t{j} = nan(n_MC, 1);
end
use_fault_depth = strcmp(G.fault_stress_mode, 'fault_depth') && ...
    isfield(FAULTS, 'fault_coord_z') && ...
    ~isempty(FAULTS.fault_coord_z);
for i = 1:n_MC
    if stress_uncertain == "yes"
        SHdir = ss.r_dir(i);
        SHmax = ss.r_SHmax(i);
        Shmin = ss.r_Shmin(i);
        Sv    = ss.r_Sv(i);
        p0    = ss.r_p(i);
    else
        SHdir = mean_stress.SHdir;
        SHmax = mean_stress.SHmax;
        Shmin = mean_stress.Shmin;
        Sv    = mean_stress.Sv;
        p0    = mean_stress.p;
    end

    if fault_uncertainty == "yes"
        mu_i = r_mu(i);
    else
        mu_i = G.fault_friction;
    end

    for j = 1:nr_fault
        if fault_uncertainty == "yes"
            dip_ij = r_dip(i,j);
            azi_ij = r_azi(i,j);
        else
            dip_ij = FAULTS.fault_dip(j);
            azi_ij = FAULTS.fault_azi(j);
        end
        if use_fault_depth
            stress_base = fault_stress_state(G, FAULTS.fault_coord_z(j));
            SHmax_j = stress_base.SHmax + (SHmax - mean_stress.SHmax);
            Shmin_j = stress_base.Shmin + (Shmin - mean_stress.Shmin);
            Sv_j    = stress_base.Sv    + (Sv    - mean_stress.Sv);
            p0_j    = stress_base.p     + (p0    - mean_stress.p);
        else
            SHmax_j = SHmax;
            Shmin_j = Shmin;
            Sv_j    = Sv;
            p0_j    = p0;
        end
        [Sn, Tau] = stress_projection(SHdir, SHmax_j, Shmin_j, Sv_j, azi_ij, dip_ij);
        deltap_cr_s{j}(i) = Sn - p0_j  - Tau / mu_i;
        deltap_cr_t{j}(i) = Sn - p0_j  + tens;
    end
end

samples.deltap_cr_s = deltap_cr_s;
samples.deltap_cr_t = deltap_cr_t;
end


function GEO = aggregate_raw_samples(GEO, samples, FAULTS)
nr_fault = FAULTS.nr_fault;
deltap_cr_samples = cell(1, nr_fault);

for j = 1:nr_fault
    vals = min(samples.deltap_cr_s{j}, samples.deltap_cr_t{j});
    vals(vals <= 0) = NaN;
    deltap_cr_samples{j} = vals;
end

GEO.deltap_cr_samples = deltap_cr_samples;
end


function FAULTS = get_fault_geometry_deterministic(G)
FAULTS.nr_fault      = G.nr_fault;
FAULTS.fault_coord_x = G.fault_coord_x;
FAULTS.fault_coord_y = G.fault_coord_y;
FAULTS.fault_dip     = G.fault_dip;
FAULTS.fault_azi     = G.fault_azi;
if isfield(G, 'fault_length')
    FAULTS.fault_length = G.fault_length;
end
if isfield(G, 'fault_coord_z') && ~isempty(G.fault_coord_z)
    FAULTS.fault_coord_z = G.fault_coord_z;
end
end


function FAULTS = get_fault_geometry_stochastic(G, area_res, x_length, y_length, thickness)
V     = area_res * 1e6 * thickness;
P32   = 0.0001;
alpha = 3;
Lmin  = 800;
Lmax  = 4000;

if alpha ~= 1
    C = (1 - alpha) / (Lmax^(1-alpha) - Lmin^(1-alpha));
else
    C = 1 / log(Lmax / Lmin);
end
if alpha ~= 3
    mean_A = C * (pi/4) * (Lmax^(3-alpha) - Lmin^(3-alpha)) / (3 - alpha);
else
    mean_A = C * (pi/4) * log(Lmax / Lmin);
end

Nv       = P32 / mean_A;
nr_fault = round(Nv * V);

U = rand(nr_fault, 1);
fault_length = (Lmin^(1-alpha) + U .* (Lmax^(1-alpha) - Lmin^(1-alpha))).^(1/(1-alpha));

FAULTS.nr_fault      = nr_fault;
FAULTS.fault_coord_x = rand(nr_fault, 1) * x_length;
FAULTS.fault_coord_y = rand(nr_fault, 1) * y_length;
FAULTS.fault_length  = fault_length;
FAULTS.fault_dip     = 60 * ones(nr_fault, 1);
FAULTS.fault_azi     = 45 * ones(nr_fault, 1);
end


function pd = build_dist(dist_type, offset, depth, grad, a_grad) %#ok<DEFNU>
mean_val = offset + depth * grad;
half     = depth * a_grad;
if strcmp(dist_type, 'Uni')
    pd = makedist('Uniform', 'lower', mean_val-half, 'upper', mean_val+half);
else
    pd = makedist('Normal', 'mu', mean_val, 'sigma', half);
end
end

function pd = build_dist_angle(dist_type, center, half_range)
if strcmp(dist_type, 'Uni')
    pd = makedist('Uniform', 'lower', center-half_range, 'upper', center+half_range);
else
    pd = makedist('Normal', 'mu', center, 'sigma', half_range);
end
end