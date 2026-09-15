function [P, GEOM] = compute_pressure_buildup(S, R, GEO)
% Well grid generation and superposed overpressure at centre/fault locations.
% GEOM caches perm-independent geometry for reuse in MC and plots.

% local shortcuts used repeatedly in the Nordbotten loop
perm    = R.perm;
thick   = R.thick;
por     = R.por;
visc_w  = R.visc_w;
visc_c  = R.visc_c;
compr   = R.compr;
rc      = R.rc;
gamma   = R.gamma;
delta   = R.delta;
omega   = R.omega;
use_sf  = strcmp(R.use_shape_factor, 'yes');
x_length = R.x_length;
y_length = R.y_length;

correction        = S.correction;
center_point_type = S.center_point_type;
rw                = S.rw;
time_yr           = S.time_yr;
nr_dist           = S.nr_dist;
nr_well_max       = S.nr_well_max;
dist_min          = S.dist_min;
dist_max          = S.dist_max;

nr_fault      = GEO.nr_fault;
fault_coord_x = GEO.fault_coord_x;
fault_coord_y = GEO.fault_coord_y;

% check all faults lie within the closed reservoir domain
no_faults_in_domain = false;
if center_point_type == "Fault" && nr_fault > 0 && isfinite(rc)
    tol = 0.01 * min(x_length, y_length);   % 1% tolerance for edge faults
    inside = (fault_coord_x >= -tol) & (fault_coord_x <= x_length + tol) & ...
        (fault_coord_y >= -tol) & (fault_coord_y <= y_length + tol);
    if ~any(inside)
        no_faults_in_domain = true;
        fprintf(['\n  *** WARNING: No faults detected within the closed reservoir domain.\n' ...
            '      Reservoir bounds: [0, %.0f] x [0, %.0f] m\n' ...
            '      All %d faults lie outside — fault reactivation constraint removed.\n' ...
            '      Falling back to central well as the reference point for capacity.\n\n'], ...
            x_length, y_length, nr_fault);
    end
end

time        = time_yr * 86400 * 365;
R_influence = sqrt(2.246 * perm * time / (visc_w * compr));  % Theis radius of influence
M0          = perm / 1e-13;

if use_sf
    fprintf('compute_pressure_buildup: shape factor ON | X_e=%.0fm Y_e=%.0fm rc=%.0fm\n', ...
        x_length, y_length, rc);
end

% resolve auto nr_well_max and dist_max
nr_well_max_is_auto = ischar(nr_well_max) && strcmp(nr_well_max,'auto');
dist_max_is_auto    = ischar(dist_max)    && strcmp(dist_max,   'auto');

if nr_well_max_is_auto
    nx_max = ceil(x_length / (dist_min*1000));
    ny_max = ceil(y_length / (dist_min*1000));
    nr_well_max = nx_max * ny_max;
end
if dist_max_is_auto
    dist_max = 0.5 * sqrt(x_length^2 + y_length^2) / 1000;   % half-diagonal [km]
end

% input checks
if ~nr_well_max_is_auto && nr_well_max < 1
    error('CO2BLOCK:invalidNrWellMax', ...
        ['nr_well_max is set to %d in General_settings.xlsx, but it must be >= 1.\n' ...
        'Please enter a positive integer (or ''auto'') for nr_well_max.'], ...
        nr_well_max);
end
if nr_well_max_is_auto && nr_well_max < 1
    error('CO2BLOCK:invalidDistMin', ...
        ['dist_min (%.4g km) is too large for the reservoir dimensions ' ...
        '(%.0f m x %.0f m).\n\n' ...
        'With nr_well_max=auto, at least one well requires:\n' ...
        '  dist_min <= min(x_length, y_length)/1000 = %.4g km\n\n' ...
        'Please reduce dist_min in General_settings.xlsx,\n' ...
        'or increase the reservoir dimensions.'], ...
        dist_min, x_length, y_length, min(x_length, y_length)/1000);
end
if nr_dist < 1
    error('CO2BLOCK:invalidNrDist', ...
        ['nr_dist is set to %d in General_settings.xlsx, but it must be >= 1.\n' ...
        'Please enter a positive integer for nr_dist.'], ...
        nr_dist);
end
if nr_dist < 2
    warning('CO2BLOCK:fewDistanceSteps', ...
        ['nr_dist=%d will produce only 1 distance column.\n' ...
        'Contour plots require at least a 2x2 data matrix and will\n' ...
        'be skipped. Set nr_dist >= 2 for full contour output.'], ...
        nr_dist);
end
if ~nr_well_max_is_auto && nr_well_max <= 2
    warning('CO2BLOCK:fewWellScenarios', ...
        ['nr_well_max=%d will produce only 1 well scenario.\n' ...
        'Contour plots require at least a 2x2 data matrix and will\n' ...
        'be skipped. Set nr_well_max >= 3 (or ''auto'') for full contour output.'], ...
        nr_well_max);
end
if dist_min >= dist_max
    if dist_max_is_auto
        dist_max_auto = 0.5 * sqrt(x_length^2 + y_length^2) / 1000;
        error('CO2BLOCK:invalidDistMin', ...
            ['dist_min (%.4g km) is too large for the reservoir ' ...
            '(%.0f m x %.0f m).\n\n' ...
            'With dist_max=auto, the maximum inter-well distance is\n' ...
            'the half-diagonal of the rectangle:\n' ...
            '  dist_max = 0.5*sqrt(%.0f^2 + %.0f^2)/1000 = %.4g km\n\n' ...
            'dist_min must be strictly less than dist_max.\n' ...
            'Please set dist_min < %.4g km in General_settings.xlsx,\n' ...
            'or increase the reservoir dimensions.'], ...
            dist_min, x_length, y_length, ...
            x_length, y_length, dist_max_auto, dist_max_auto);
    else
        error('CO2BLOCK:invalidDistMin', ...
            ['dist_min (%.4g km) must be strictly less than dist_max (%.4g km).\n\n' ...
            'Please fix dist_min or dist_max in General_settings.xlsx.'], ...
            dist_min, dist_max);
    end
end

d_list = linspace(dist_min, dist_max, nr_dist);

if center_point_type == "Desired"
    center_coord_x = S.desired_x * 1000;
    center_coord_y = S.desired_y * 1000;
end

% preallocate
b                = [];
well_list        = [];
d_max_vec        = [];
x_grid_list      = [];
y_grid_list      = [];
distance_matrix  = [];
p_sup_vec        = [];
Q0_vec           = [];
well_coords_x    = {};
well_coords_y    = {};
dist_vec_cache   = {};
dist_fault_cache = {};
shape_CA_cache   = {};
shape_Rext_cache = {};
p_fault          = {};
p_c_last         = NaN;
w_id             = 0;

% generate valid (nx, ny) well grid configurations
grid_pairs  = [];
nx_max      = ceil(x_length / (dist_min*1000));
ny_max      = ceil(y_length / (dist_min*1000));
x_is_longer = x_length > y_length;
for nx = 1:nx_max
    for ny = 1:ny_max
        if nx * ny > nr_well_max, continue; end
        if x_is_longer && nx < ny, continue; end
        if ~x_is_longer && ny < nx, continue; end
        grid_pairs = [grid_pairs; nx, ny, nx*ny]; %#ok<AGROW>
    end
end
grid_pairs = sortrows(grid_pairs, [3 1 2]);

% main loop over grid configurations and well spacings
for gp = 1:size(grid_pairs, 1)
    x_grid_num = grid_pairs(gp, 1);
    y_grid_num = grid_pairs(gp, 2);

    w_id = w_id + 1;
    w    = x_grid_num * y_grid_num;
    well_list(w_id)   = w;
    x_grid_list(w_id) = x_grid_num;
    y_grid_list(w_id) = y_grid_num;

    d_max_vec(w_id) = sqrt(x_length * y_length / w) / 1000;   % sqrt(area per well) [km]

    Q0           = M0 * 1e9 / R.dens_c / 365 / 86400 / w;
    Q0_vec(w_id) = Q0;
    csi          = sqrt(Q0 * time / pi / por / thick);
    psi          = exp(omega) * csi;
    p_c_last     = (Q0 * visc_w) / (2*pi*thick*perm) / 1e6;
    b(w_id)      = (visc_w - visc_c) / (4*pi*perm*thick);

    for d = 1:nr_dist
        distance_km = d_list(d);
        distance    = distance_km * 1000;

        if distance_km > d_max_vec(w_id)
            distance_matrix(w_id, d) = NaN;
            p_sup_vec(w_id, d)       = NaN;
            dist_vec_cache{w_id, d}  = [];
            dist_fault_cache{w_id, d}= [];
            continue;
        end
        distance_matrix(w_id, d) = distance;

        offset_x = (x_length - (x_grid_num-1)*distance) / 2;
        offset_y = (y_length - (y_grid_num-1)*distance) / 2;

        if offset_x < 0 || offset_y < 0
            distance_matrix(w_id, d) = NaN;
            p_sup_vec(w_id, d)       = NaN;
            dist_vec_cache{w_id, d}  = [];
            dist_fault_cache{w_id, d}= [];
            continue;
        end

        margin = min(offset_x, offset_y);
        if x_grid_num > 1
            dx = (x_length - 2*margin) / (x_grid_num - 1);
            x_start = margin;
        else
            dx = 0;
            x_start = x_length / 2;
        end
        if y_grid_num > 1
            dy = (y_length - 2*margin) / (y_grid_num - 1);
            y_start = margin;
        else
            dy = 0;
            y_start = y_length / 2;
        end
        xv = x_start + (0:x_grid_num-1) * dx;
        yv = y_start + (y_grid_num-1:-1:0) * dy;
        [wells_x, wells_y] = meshgrid(xv, yv);
        well_coords_x{w_id, d} = wells_x;
        well_coords_y{w_id, d} = wells_y;

        if use_sf
            [C_A_vec, R_ext_vec] = compute_shape_factors( ...
                wells_x(:), wells_y(:), x_length, y_length);
            shape_CA_cache{w_id, d}   = C_A_vec;
            shape_Rext_cache{w_id, d} = R_ext_vec;
        else
            shape_CA_cache{w_id, d}   = [];
            shape_Rext_cache{w_id, d} = [];
        end

        if center_point_type == "Central" || no_faults_in_domain
            center_coord_x = x_length / 2;
            center_coord_y = y_length / 2;
        end

        if center_point_type == "Central" || center_point_type == "Desired" ...
                || no_faults_in_domain
            dv = sqrt((wells_x - center_coord_x).^2 + ...
                (wells_y - center_coord_y).^2);
            dv(dv == 0) = rw;
            dist_vec_cache{w_id, d} = dv(:);
        else
            dist_vec_cache{w_id, d} = [];
        end

        if center_point_type == "Fault" && nr_fault > 0 && ~no_faults_in_domain
            wx = wells_x(:);
            wy = wells_y(:);
            nw = numel(wx);
            df = zeros(nr_fault, nw);
            for f = 1:nr_fault
                r_f = sqrt((wx - fault_coord_x(f)).^2 + ...
                    (wy - fault_coord_y(f)).^2);
                r_f(r_f == 0) = rw;
                df(f,:) = r_f;
            end
            dist_fault_cache{w_id, d} = df;
        else
            dist_fault_cache{w_id, d} = [];
        end

        apply_corr = strcmp(correction,'on') && w >= 9 && ...
            R_influence * csi / distance^2 >= 1;
        if strcmp(correction,'on') && apply_corr
            b(w_id) = (visc_w - visc_c) / (4*pi*perm*thick*(1+w/4));
        end
        sup_error = 0;
        if apply_corr
            sup_error = w * delta / 4 * log(R_influence * csi / distance^2);
        end

        dv = dist_vec_cache{w_id, d};
        if ~isempty(dv)
            PD = Nordbotten_solution(dv, R_influence, psi, rc, gamma, ...
                use_sf, shape_CA_cache{w_id,d}, shape_Rext_cache{w_id,d}, R.area_res);
            p_sup_vec(w_id, d) = sum(PD)*p_c_last - sup_error*p_c_last;
        end

        if center_point_type == "Fault" && nr_fault > 0 && ~no_faults_in_domain
            df     = dist_fault_cache{w_id, d};
            pf_vec = zeros(nr_fault, 1);
            C_A_v  = shape_CA_cache{w_id, d};
            Rx_v   = shape_Rext_cache{w_id, d};
            for f = 1:nr_fault
                PD_f    = Nordbotten_solution(df(f,:)', R_influence, psi, rc, gamma, ...
                    use_sf, C_A_v, Rx_v, R.area_res);
                p_local = sum(PD_f) * p_c_last;
                if apply_corr, p_local = p_local - sup_error*p_c_last; end
                pf_vec(f) = p_local;
            end
            if isfinite(rc)
                tol = 0.01 * min(x_length, y_length);
                outside = (fault_coord_x < -tol) | (fault_coord_x > x_length + tol) | ...
                    (fault_coord_y < -tol) | (fault_coord_y > y_length + tol);
                pf_vec(outside) = 0;
            end
            p_fault{w_id}{d} = pf_vec;
        end

    end % distance loop
end % grid_pairs loop

% geometry cache — perm-independent, reused by MC loop and plots
GEOM.well_coords_x = well_coords_x;
GEOM.well_coords_y = well_coords_y;
GEOM.dist_vec_cache = dist_vec_cache;
GEOM.dist_fault_cache = dist_fault_cache;
GEOM.shape_CA_cache = shape_CA_cache;
GEOM.shape_Rext_cache = shape_Rext_cache;
GEOM.use_sf = use_sf;
GEOM.x_length = x_length;
GEOM.y_length = y_length;
GEOM.d_list = d_list;
GEOM.well_list = well_list;
GEOM.d_max = d_max_vec;
GEOM.distance_matrix = distance_matrix;
GEOM.x_grid_list = x_grid_list;
GEOM.y_grid_list = y_grid_list;
GEOM.time = time;
GEOM.nr_fault = nr_fault;
GEOM.rw = rw;
GEOM.correction = correction;
GEOM.center_point_type = center_point_type;
GEOM.R_influence = R_influence;
GEOM.no_faults_in_domain = no_faults_in_domain;

% pressure output
P.d_list = d_list;
P.well_list = well_list;
P.d_max = d_max_vec;
P.x_grid_list = x_grid_list;
P.y_grid_list = y_grid_list;
P.distance_matrix = distance_matrix;
P.x_length = x_length;
P.y_length = y_length;
P.well_coords_x = well_coords_x;
P.well_coords_y = well_coords_y;
P.b = b;
P.Q0_vec = Q0_vec;
P.p_c = p_c_last;
P.M0 = M0;
P.p_sup_vec = p_sup_vec;
P.p_fault = p_fault;
P.p_sup_2Dgrid = {};   % filled on demand in CO2BLOCK_plots
P.Mesh_grid = {};
P.no_faults_in_domain = no_faults_in_domain;

fprintf('compute_pressure_buildup: %d well scenarios, %d distances\n', ...
    length(well_list), nr_dist);
end


function [C_A, R_ext_well] = compute_shape_factors(wx, wy, X_e, Y_e)
% Dietz shape factor C_A and max corner distance R_ext for each well.
f1 = X_e / Y_e;
f2 = wx / X_e;
f3 = wy / Y_e;
f4 = (8*pi*f1) .* (1/6 - f2/2 + f2.^2/2);
C_A = 88.6657 * f1 * sin(pi*f3).^2 ./ exp(f4);

persistent sf_printed;
if isempty(sf_printed)
    fprintf('  compute_shape_factors: nw=%d | f2=[%.3f..%.3f] f3=[%.3f..%.3f] C_A=[%.3f..%.3f]\n', ...
        numel(wx), min(f2), max(f2), min(f3), max(f3), min(C_A), max(C_A));
    sf_printed = true;
end

corners = [0, 0; X_e, 0; 0, Y_e; X_e, Y_e];
nw      = numel(wx);
R_ext_well = zeros(nw, 1);
for iw = 1:nw
    dx = corners(:,1) - wx(iw);
    dy = corners(:,2) - wy(iw);
    R_ext_well(iw) = max(sqrt(dx.^2 + dy.^2));
end
end