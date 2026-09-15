function CR = compute_reactivation_risk(S, P, GEO, MC)
% Criticality ratio Dp_inj / Dp_cr for each fault.
% Deterministic pass at reference perm + optional per-realization pass for MC.
% Tie-breaking when faults share the same max ratio is deferred to
% compute_storage_capacity, which needs Q to resolve it anyway.

if nargin < 4, MC = []; end

center_point_type = S.center_point_type;
nr_fault      = GEO.nr_fault;
deltap_cr_det = GEO.deltap_cr_det;
well_list = P.well_list;
d_list    = P.d_list;
nr_dist   = numel(d_list);
p_fault   = P.p_fault;
nw        = numel(well_list);

have_dc_samples = isfield(GEO,'deltap_cr_samples') && ...
                  iscell(GEO.deltap_cr_samples)    && ...
                  numel(GEO.deltap_cr_samples) == nr_fault;

% early exit — non-fault mode or no faults present
if center_point_type ~= "Fault" || nr_fault == 0
    CR.crit_ratio_det            = {};
    CR.max_ratio_det             = [];
    CR.p1_fault_det              = [];
    CR.p2_fault_det              = [];
    CR.tied_det                  = {};
    CR.most_critical_fault_table = table();
    CR.max_ratio_P10 = [];  CR.max_ratio_P50 = [];  CR.max_ratio_P90 = [];
    CR.p1_samp_raw = [];    CR.p2_samp_raw = [];    CR.mr_samp_raw = [];
    CR.fstar_samp  = [];    CR.tied_samp   = {};
    return;
end

% --- deterministic pass (reference permeability, mean stress) ---

crit_ratio_det = {};
max_ratio_det  = NaN(nw, nr_dist);
p1_fault_det   = NaN(nw, nr_dist);
p2_fault_det   = NaN(nw, nr_dist);
tied_det       = cell(nw, nr_dist);

tbl = table([],[],[],[],[],[],[],[], 'VariableNames', ...
    {'w_id','wells','d_km','coord_x','coord_y','max_crit_ratio','delta_p','deltap_cr'});

for wi = 1:nw
    if isempty(p_fault) || wi > numel(p_fault) || isempty(p_fault{wi}), continue; end
    for di = 1:nr_dist
        if di > numel(p_fault{wi}) || isempty(p_fault{wi}{di}), continue; end
        pf = p_fault{wi}{di};

        r = pf ./ deltap_cr_det(:);
        r(deltap_cr_det(:) <= 0) = NaN;
        crit_ratio_det{wi}{di} = r;

        mr = max(r);
        if ~isfinite(mr), continue; end

        % collect all faults at the maximum — tie-breaking happens later
        tied = find(abs(r - mr) < 1e-12 & isfinite(r));

        max_ratio_det(wi,di) = mr;
        p1_fault_det(wi,di)  = pf(tied(1)) * 1e6;
        p2_fault_det(wi,di)  = deltap_cr_det(tied(1)) * 1e6;
        tied_det{wi,di}      = tied;

        tbl = [tbl; {wi, well_list(wi), d_list(di), ...
            GEO.fault_coord_x(tied(1)), GEO.fault_coord_y(tied(1)), ...
            mr, pf(tied(1)), deltap_cr_det(tied(1))}]; %#ok<AGROW>
    end
end

CR.crit_ratio_det            = crit_ratio_det;
CR.max_ratio_det             = max_ratio_det;
CR.p1_fault_det              = p1_fault_det;
CR.p2_fault_det              = p2_fault_det;
CR.tied_det                  = tied_det;
CR.most_critical_fault_table = tbl;

% --- MC pass (skipped when neither perm nor Dp_cr uncertainty is active) ---

have_perm_mc = ~isempty(MC) && isfield(MC,'perm_samples') && ~isempty(MC.perm_samples);

if ~have_perm_mc && ~have_dc_samples
    CR.max_ratio_P10 = [];  CR.max_ratio_P50 = [];  CR.max_ratio_P90 = [];
    CR.p1_samp_raw = [];    CR.p2_samp_raw = [];    CR.mr_samp_raw = [];
    CR.fstar_samp  = [];    CR.tied_samp   = {};
    return;
end

if have_perm_mc
    n_MC  = numel(MC.perm_samples);
    k_ref = P.perm_ref;
else
    n_MC  = numel(GEO.deltap_cr_samples{1});
    k_ref = [];
end

mr_samp   = NaN(n_MC, nw, nr_dist);
p1_samp   = NaN(n_MC, nw, nr_dist);
p2_samp   = NaN(n_MC, nw, nr_dist);
fs_samp   = NaN(n_MC, nw, nr_dist);
tied_samp = cell(n_MC, nw, nr_dist);

for i = 1:n_MC

    pf_scale = 1.0;
    if have_perm_mc
        pf_scale = k_ref / MC.perm_samples(i);
    end

    if have_dc_samples
        dc_i = cellfun(@(s) s(i), GEO.deltap_cr_samples);
        dc_i = dc_i(:);
    else
        dc_i = deltap_cr_det(:);
    end

    for wi = 1:nw
        if isempty(p_fault) || wi > numel(p_fault) || isempty(p_fault{wi}), continue; end
        for di = 1:nr_dist
            if di > numel(p_fault{wi}) || isempty(p_fault{wi}{di}), continue; end

            pf_i = p_fault{wi}{di} * pf_scale;
            r_i  = pf_i ./ dc_i;
            r_i(dc_i <= 0 | ~isfinite(dc_i)) = NaN;

            mr = max(r_i);
            if ~isfinite(mr), continue; end

            tied = find(abs(r_i - mr) < 1e-12 & isfinite(r_i));

            mr_samp(i,wi,di)   = mr;
            p1_samp(i,wi,di)   = p_fault{wi}{di}(tied(1));
            p2_samp(i,wi,di)   = dc_i(tied(1));
            fs_samp(i,wi,di)   = tied(1);
            tied_samp{i,wi,di} = tied;
        end
    end
end

CR.max_ratio_P10 = squeeze(prctile(mr_samp, 10, 1));
CR.max_ratio_P50 = squeeze(prctile(mr_samp, 50, 1));
CR.max_ratio_P90 = squeeze(prctile(mr_samp, 90, 1));
CR.p1_samp_raw   = p1_samp;
CR.p2_samp_raw   = p2_samp;
CR.mr_samp_raw   = mr_samp;
CR.fstar_samp    = fs_samp;
CR.tied_samp     = tied_samp;

fprintf('  reactivation risk: %d realizations done\n', n_MC);
end