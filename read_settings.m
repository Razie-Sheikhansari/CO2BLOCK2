function S = read_settings(fname)
% Read run-control parameters from General_settings.xlsx.

    if nargin < 1
        fname = 'General_settings.xlsx';
    end

    %% read file into parameter map
    raw    = readcell(fname, 'Sheet', 'General_settings', 'Range', 'A:B');
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

    %% well configuration and injection parameters
    S.rw          = num('rw');
    S.dist_min    = num('dist_min');
    S.nr_dist     = num('nr_dist');
    S.time_yr     = num('time_yr');
    S.maxQ        = num('maxQ');
    S.correction  = str('correction');
    S.nr_well_max = parse_auto(params('nr_well_max'));
    S.dist_max    = parse_auto(params('dist_max'));

    %% well placement reference point
    S.center_point_type = string(str('Reference Point'));
    S.desired_x         = num('Desired x');
    S.desired_y         = num('Desired y');

    % fault_type is only required when center_point_type = "Fault"
    if S.center_point_type == "Fault"
        S.fault_type = string(str('Fault Type'));
    else
        if has('Fault Type') && ~cell_is_blank(params('Fault Type'))
            S.fault_type = string(str('Fault Type'));
        else
            S.fault_type = "N/A";
        end
    end

    %% stress_uncertain
    if has('stress_uncertain') && ~cell_is_blank(params('stress_uncertain'))
        S.stress_uncertain = string(str('stress_uncertain'));
    else
        S.stress_uncertain = "no";
    end

    %% fault_uncertainty
    if has('fault_uncertainty') && ~cell_is_blank(params('fault_uncertainty'))
        S.fault_uncertainty = string(str('fault_uncertainty'));
    else
        S.fault_uncertainty = "no";
    end

    %% perm_uncertain
    if has('perm_uncertain') && ~cell_is_blank(params('perm_uncertain'))
        S.perm_uncertain = string(str('perm_uncertain'));
    else
        S.perm_uncertain = "no";
    end

    %% n_MC
    if has('n_MC') && ~cell_is_blank(params('n_MC'))
        S.n_MC     = num('n_MC');
        S.n_MC_set = true;
    else
        S.n_MC     = NaN;
        S.n_MC_set = false;
    end

    %% validate
    if ~ismember(S.stress_uncertain, ["yes","no"])
        error('read_settings: stress_uncertain must be "yes" or "no", got "%s".', ...
            S.stress_uncertain);
    end
    if ~ismember(S.fault_uncertainty, ["yes","no"])
        error('read_settings: fault_uncertainty must be "yes" or "no", got "%s".', ...
            S.fault_uncertainty);
    end
    if ~ismember(S.perm_uncertain, ["yes","no"])
        error('read_settings: perm_uncertain must be "yes" or "no", got "%s".', ...
            S.perm_uncertain);
    end
    if S.n_MC_set && S.n_MC < 1
        error('read_settings: n_MC must be >= 1, got %d.', S.n_MC);
    end

    any_uncertain = S.stress_uncertain == "yes" || ...
                    S.fault_uncertainty == "yes" || ...
                    S.perm_uncertain    == "yes";
    if any_uncertain && ~S.n_MC_set
        error(['read_settings: n_MC is blank in General_settings.xlsx,\n' ...
               'but at least one uncertainty source is active.\n' ...
               'Please set n_MC to a positive integer (e.g. 1000).']);
    end

    if S.fault_uncertainty == "yes" && S.center_point_type ~= "Fault"
        warning('read_settings: fault_uncertainty=yes is ignored when center_point_type="%s".', ...
            S.center_point_type);
    end
end


function out = parse_auto(val)
    if ischar(val) && strcmpi(strtrim(val), 'auto')
        out = 'auto';
    else
        out = double(val);
    end
end


function tf = cell_is_blank(val)
    if ismissing(val)
        tf = true;
    elseif isnumeric(val) && isnan(val)
        tf = true;
    else
        s  = strtrim(char(string(val)));
        tf = isempty(s) || strcmp(s, '-');
    end
end