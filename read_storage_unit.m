function R = read_storage_unit(path, name)
% Read reservoir/storage unit properties from Storage_unit_properties.xlsx.
% Does not call eos; fluid properties are computed later by compute_fluid_properties.

    if nargin < 1, path = ''; end
    if nargin < 2, name = 'Storage_unit_properties.xlsx'; end

    %% read file into parameter map
    raw    = readcell(fullfile(path, name), 'Sheet', 1, 'Range', 'A:B');
    raw    = raw(cellfun(@(x) ischar(x) && ~isempty(strtrim(x)), raw(:,1)), :);
    params = containers.Map(strtrim(raw(:,1)), raw(:,2));

    function v = num(key)
        v = double(params(key));
    end
    function v = str(key)
        v = strtrim(char(params(key)));
    end
    function v = has(key)
        v = isKey(params, key);
    end

    %% reservoir geometry
    domain_type  = str('domain_type');
    thick        = num('thick');
    area_res     = num('area_res');

    %% flow properties
    perm         = num('perm') * 1e-15;      % mD -> m2
    por          = num('por');
    cr           = num('cr')  / 1e6;         % MPa-1 -> Pa-1  (0 -> default in compute_fluid_properties)
    cw           = num('cw')  / 1e6;         % MPa-1 -> Pa-1  (0 -> default in compute_fluid_properties)

    %% fluid identity (raw; eos computed later by compute_fluid_properties)
    dens_c       = num('dens_c') * 1e3;      % ton/m3 -> kg/m3  (0 = compute from eos)
    visc_c       = num('visc_c') / 1e3;      % cp -> Pa.s        (0 = compute from eos)
    visc_w       = num('visc_w') / 1e3;      % cp -> Pa.s        (0 = compute from eos)
    salinity     = num('salinity') / 1e6;    % ppm -> fraction

    %% salinity default (no P or T needed)
    if salinity == 0 || isnan(salinity)
        salinity = 180000 / 1e6;             % 180000 ppm -> fraction
    end

    %% domain type to outer radius
    switch lower(domain_type)
        case 'open';   rc = inf;
        case 'closed'; rc = sqrt(area_res * 1e6 / pi);
        otherwise;     error('read_storage_unit: unknown domain_type "%s"', domain_type);
    end

    %% pack reservoir fields
    % compr, gamma, delta, omega are not set here; computed by compute_fluid_properties.
    R.domain_type = lower(domain_type);   % 'open' or 'closed'
    R.thick    = thick;     R.area_res = area_res;  R.rc       = rc;
    R.perm     = perm;      R.por      = por;
    R.cr       = cr;        R.cw       = cw;
    R.dens_c   = dens_c;    R.visc_c   = visc_c;    R.visc_w   = visc_w;
    R.salinity = salinity;

    %% permeability uncertainty distribution
    % perm_uncertain and n_MC are read from General_settings; only shape parameters remain here.
    if has('perm_dist')
        R.perm_dist = str('perm_dist');             % 'Nor' or 'Uni'
    else
        R.perm_dist = 'Nor';
    end
    if has('perm_std')
        R.perm_std = num('perm_std');
    else
        R.perm_std = 0.5;                           % default +/-0.5 log10 units
    end

    %% shape factor
    % Optional rows: use_shape_factor, x_length, y_length. Default: no.
    if has('use_shape_factor')
        R.use_shape_factor = lower(str('use_shape_factor'));  % 'yes' or 'no'
    else
        R.use_shape_factor = 'no';
    end

    if strcmp(R.use_shape_factor, 'yes')
        if ~has('x_length') || ~has('y_length')
            error('read_storage_unit: use_shape_factor=yes requires x_length and y_length rows in Excel.');
        end
        R.x_length = num('x_length') * 1000;   % km -> m
        R.y_length = num('y_length') * 1000;   % km -> m
        fprintf('read_storage_unit: use_shape_factor=YES | x_length=%.4g km  y_length=%.4g km\n', ...
            R.x_length/1000, R.y_length/1000);
    else
        if has('x_length') && has('y_length')
            R.x_length = num('x_length') * 1000;   % km -> m
            R.y_length = num('y_length') * 1000;   % km -> m
            fprintf('read_storage_unit: use_shape_factor=NO | x_length=%.4g km  y_length=%.4g km\n', ...
                R.x_length/1000, R.y_length/1000);
        else
            R.x_length = sqrt(R.area_res) * 1000;  % fallback square [m]
            R.y_length = sqrt(R.area_res) * 1000;
            fprintf('read_storage_unit: use_shape_factor=NO | fallback square: x_length=%.4g km  y_length=%.4g km\n', ...
                R.x_length/1000, R.y_length/1000);
        end
    end
end