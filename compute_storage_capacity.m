function OUT = compute_storage_capacity(S, R, P, CR, MC, GEO, GEOM)
% Storage capacity: deterministic, fault-constrained, and joint MC percentiles.
% lambertw calls wrapped in real() to handle rare complex output near the branch cut at -1/e.

    time = S.time_yr * 86400 * 365;
    is_fault   = (S.center_point_type == "Fault") && ...
                 ~(isfield(GEOM,'no_faults_in_domain') && GEOM.no_faults_in_domain);
    is_desired = (S.center_point_type == "Desired");

    % deterministic limiting Dp_cr
    if is_fault
        dc_det = min(GEO.deltap_cr_det(GEO.deltap_cr_det > 0), [], 'omitnan');
    else
        dc_det = GEO.deltap_cr_det;
        if ~isscalar(dc_det), dc_det = min(dc_det(dc_det > 0)); end
    end

    % resolve tied faults before capacity solve (deterministic)
    if is_fault && isfield(CR,'tied_det') && ~isempty(CR.tied_det)
        for wi = 1:numel(P.well_list)
            for di = 1:numel(P.d_list)
                if isempty(CR.tied_det{wi,di}) || numel(CR.tied_det{wi,di}) <= 1
                    continue;
                end
                tied  = CR.tied_det{wi,di};
                b_wi  = P.b(wi);
                Q0_wi = P.Q0_vec(wi);
                q_tied = NaN(numel(tied), 1);
                for t = 1:numel(tied)
                    f   = tied(t);
                    pf1 = P.p_fault{wi}{di}(f) * 1e6;
                    pf2 = GEO.deltap_cr_det(f) * 1e6;
                    lw  = -pf2 / Q0_wi / b_wi * exp(-pf1 / Q0_wi / b_wi);
                    if lw >= -1/exp(1) && lw < 0 && isfinite(lw)
                        q_tied(t) = real(-pf2 / b_wi / lambertw(-1, lw));
                    end
                end
                [~, best] = min(q_tied);
                if ~isnan(q_tied(best))
                    winner = tied(best);
                    CR.p1_fault_det(wi,di) = P.p_fault{wi}{di}(winner) * 1e6;
                    CR.p2_fault_det(wi,di) = GEO.deltap_cr_det(winner) * 1e6;
                end
            end
        end
    end

    [Q_M_each_det, V_M_det, Table_Q_det, Table_V_det, ...
     Q_M_each_fault_det, V_M_fault_det, Table_Q_fault_det, Table_V_fault_det] = ...
        capacity_one_mode(P.d_list, P.well_list, P.x_grid_list, P.y_grid_list, ...
            P.d_max, P.p_sup_vec, P.b, P.Q0_vec, ...
            dc_det, R.por, R.thick, S.time_yr, time, R.dens_c, S.maxQ, ...
            S.center_point_type, CR.p1_fault_det, CR.p2_fault_det);

    fprintf('compute_storage_capacity: [Det] done  [%d/%d finite Q_M]\n', ...
        nnz(isfinite(Q_M_each_det(:))), numel(Q_M_each_det));
    if is_fault
        fprintf('compute_storage_capacity: [Fault-Det] done  [%d/%d finite Q_M]\n', ...
            nnz(isfinite(Q_M_each_fault_det(:))), numel(Q_M_each_fault_det));
    end

    have_perm_mc    = MC.perm_active && ~isempty(MC.perm_samples);
    have_dc_samples = isfield(GEO,'deltap_cr_samples') && ~isempty(GEO.deltap_cr_samples);
    any_mc = have_perm_mc || have_dc_samples;

    Q_M_each_joint_P10 = [];  V_M_joint_P10 = [];
    Q_M_each_joint_P50 = [];  V_M_joint_P50 = [];
    Q_M_each_joint_P90 = [];  V_M_joint_P90 = [];
    Q_M_each_joint_fault_P10 = [];  V_M_joint_fault_P10 = [];
    Q_M_each_joint_fault_P50 = [];  V_M_joint_fault_P50 = [];
    Q_M_each_joint_fault_P90 = [];  V_M_joint_fault_P90 = [];
    p1_fault_mc = [];  p2_fault_mc = [];  mc_idx_fault = [];
    Vmax_fault_samp = [];

    if any_mc
        if have_perm_mc
            n_MC = numel(MC.perm_samples);
        elseif iscell(GEO.deltap_cr_samples)
            n_MC = numel(GEO.deltap_cr_samples{1});
        else
            n_MC = numel(GEO.deltap_cr_samples);
        end

        k_ref      = R.perm;
        nw_mc      = numel(P.well_list);
        nr_dist_mc = numel(P.d_list);

        Qmax_res_samp   = NaN(n_MC, nw_mc, nr_dist_mc);
        Qmax_fault_samp = NaN(n_MC, nw_mc, nr_dist_mc);

        have_tied    = isfield(CR,'tied_samp') && ~isempty(CR.tied_samp);
        have_samp_raw = is_fault && isfield(CR,'p1_samp_raw') && ~isempty(CR.p1_samp_raw);

        for i = 1:n_MC

            if have_perm_mc
                pf_scale = k_ref / MC.perm_samples(i);
            else
                pf_scale = 1.0;
            end
            b_i  = P.b      * pf_scale;
            Q0_i = P.Q0_vec / pf_scale;

            if have_dc_samples && ~iscell(GEO.deltap_cr_samples) && numel(GEO.deltap_cr_samples) >= n_MC
                dc_i = GEO.deltap_cr_samples(i);
            elseif have_dc_samples && iscell(GEO.deltap_cr_samples)
                vals = cellfun(@(s) s(i), GEO.deltap_cr_samples);
                dc_i = min(vals(vals > 0));
                if isempty(dc_i), dc_i = NaN; end
            elseif have_dc_samples && isscalar(GEO.deltap_cr_samples)
                dc_i = GEO.deltap_cr_samples;
            else
                dc_i = dc_det;
            end

            for wi = 1:nw_mc
                b_wi  = b_i(wi);
                Q0_wi = Q0_i(wi);
                if is_fault
                    for di = 1:nr_dist_mc
                        if have_samp_raw
                            p1_wd = CR.p1_samp_raw(i, wi, di) * 1e6;
                            p2_wd = CR.p2_samp_raw(i, wi, di) * 1e6;
                        else
                            p1_wd = CR.p1_fault_det(wi, di);
                            p2_wd = dc_i * 1e6;
                        end
                        if ~isfinite(p1_wd) || ~isfinite(p2_wd) || p2_wd <= 0
                            continue;
                        end
                        lw_a = -p2_wd / Q0_wi / b_wi * exp(-p1_wd / Q0_wi / b_wi);
                        if lw_a >= -1/exp(1) && lw_a < 0 && isfinite(lw_a)
                            Qmax_res_samp(i, wi, di) = real(-p2_wd / b_wi / ...
                                lambertw(-1, lw_a)) * 86400*365*R.dens_c/1e9;
                        end
                    end
                else
                    p1_i   = P.p_sup_vec(wi, :) * 1e6;
                    p2_i   = dc_i * 1e6;
                    lw_arg = -p2_i ./ Q0_wi ./ b_wi .* exp(-p1_i ./ Q0_wi ./ b_wi);
                    ok = lw_arg >= -1/exp(1) & lw_arg < 0 & isfinite(lw_arg);
                    q_wi = NaN(1, nr_dist_mc);
                    q_wi(ok) = real(-p2_i ./ b_wi ./ lambertw(-1, lw_arg(ok)));
                    Qmax_res_samp(i, wi, :) = q_wi * 86400*365*R.dens_c/1e9;
                end
            end

            if is_fault && ~isempty(CR.fstar_samp)
                for wi = 1:nw_mc
                    if isempty(P.p_fault) || wi > length(P.p_fault), continue; end
                    b_wi  = b_i(wi);
                    Q0_wi = Q0_i(wi);
                    for di = 1:nr_dist_mc
                        if di > length(P.p_fault{wi}) || isempty(P.p_fault{wi}{di}), continue; end
                        tied = [];
                        if have_tied
                            tied = CR.tied_samp{i, wi, di};
                        end
                        if isempty(tied)
                            fs = CR.fstar_samp(i, wi, di);
                            if isnan(fs), continue; end
                            tied = fs;
                        end
                        if numel(tied) == 1
                            fs  = tied(1);
                            p1f = P.p_fault{wi}{di}(fs) * 1e6;
                            p2f = CR.p2_samp_raw(i, wi, di) * 1e6;
                            if ~isfinite(p1f) || ~isfinite(p2f) || p2f <= 0, continue; end
                            lw_arg = -p2f / Q0_wi / b_wi * exp(-p1f / Q0_wi / b_wi);
                            if lw_arg >= -1/exp(1) && lw_arg < 0 && isfinite(lw_arg)
                                q_f = real(-p2f / b_wi / lambertw(-1, lw_arg));
                                Qmax_fault_samp(i,wi,di) = q_f * 86400*365*R.dens_c/1e9;
                            end
                        else
                            q_tied = NaN(numel(tied), 1);
                            for t = 1:numel(tied)
                                f   = tied(t);
                                p1f = P.p_fault{wi}{di}(f) * 1e6;
                                if have_dc_samples && iscell(GEO.deltap_cr_samples)
                                    p2f = GEO.deltap_cr_samples{f}(i) * 1e6;
                                else
                                    p2f = GEO.deltap_cr_det(f) * 1e6;
                                end
                                if ~isfinite(p1f) || ~isfinite(p2f) || p2f <= 0, continue; end
                                lw = -p2f / Q0_wi / b_wi * exp(-p1f / Q0_wi / b_wi);
                                if lw >= -1/exp(1) && lw < 0 && isfinite(lw)
                                    q_tied(t) = real(-p2f / b_wi / lambertw(-1, lw));
                                end
                            end
                            [q_min, ~] = min(q_tied);
                            if ~isnan(q_min)
                                Qmax_fault_samp(i,wi,di) = q_min * 86400*365*R.dens_c/1e9;
                            end
                        end
                    end
                end
            end
        end  % realization loop

        % percentiles across realizations
        wm = repmat(P.well_list(:), 1, nr_dist_mc);
        Qr10 = squeeze(prctile(Qmax_res_samp, 10, 1));
        Qr50 = squeeze(prctile(Qmax_res_samp, 50, 1));
        Qr90 = squeeze(prctile(Qmax_res_samp, 90, 1));
        Q_M_each_joint_P10 = Qr10;  V_M_joint_P10 = (Qr10 .* wm) .* S.time_yr / 1000;
        Q_M_each_joint_P50 = Qr50;  V_M_joint_P50 = (Qr50 .* wm) .* S.time_yr / 1000;
        Q_M_each_joint_P90 = Qr90;  V_M_joint_P90 = (Qr90 .* wm) .* S.time_yr / 1000;

        if is_fault && any(~isnan(Qmax_fault_samp(:)))
            Qf10 = squeeze(prctile(Qmax_fault_samp, 10, 1));
            Qf50 = squeeze(prctile(Qmax_fault_samp, 50, 1));
            Qf90 = squeeze(prctile(Qmax_fault_samp, 90, 1));
            Q_M_each_joint_fault_P10 = Qf10;  V_M_joint_fault_P10 = (Qf10 .* wm) .* S.time_yr / 1000;
            Q_M_each_joint_fault_P50 = Qf50;  V_M_joint_fault_P50 = (Qf50 .* wm) .* S.time_yr / 1000;
            Q_M_each_joint_fault_P90 = Qf90;  V_M_joint_fault_P90 = (Qf90 .* wm) .* S.time_yr / 1000;
        end

        % match each percentile to the nearest MC draw so Dp_inj and Dp_cr
        % come from the same realization (prevents ratio exceeding 1)
        if is_fault && any(~isnan(Qmax_fault_samp(:)))
            pct_list     = [10 50 90];
            mc_idx_fault = NaN(nw_mc, nr_dist_mc, 3);
            p1_fault_mc  = NaN(nw_mc, nr_dist_mc, 3);
            p2_fault_mc  = NaN(nw_mc, nr_dist_mc, 3);
            for wi = 1:nw_mc
                for di = 1:nr_dist_mc
                    q  = Qmax_fault_samp(:, wi, di);
                    ok = find(isfinite(q));
                    if isempty(ok), continue; end
                    for pp = 1:3
                        target  = prctile(q(ok), pct_list(pp));
                        [~, kk] = min(abs(q(ok) - target));
                        i_star  = ok(kk);
                        mc_idx_fault(wi,di,pp) = i_star;
                        p1_fault_mc(wi,di,pp)  = CR.p1_samp_raw(i_star, wi, di) * 1e6;
                        p2_fault_mc(wi,di,pp)  = CR.p2_samp_raw(i_star, wi, di) * 1e6;
                    end
                end
            end
        end

        if is_fault && any(isfinite(Qmax_fault_samp(:)))
            wm3 = repmat(reshape(P.well_list(:), [1, nw_mc, 1]), [n_MC, 1, nr_dist_mc]);
            Vmax_fault_samp = Qmax_fault_samp .* wm3 .* S.time_yr / 1000;
        end

        fprintf('compute_storage_capacity: [Joint-MC] done (%d realizations)\n', n_MC);
    end

    % central reference baseline (always computed)
    % in Fault/Desired mode we need a fresh Central pressure field for this
    if isfield(GEO,'deltap_cr_ref') && isfinite(GEO.deltap_cr_ref) && GEO.deltap_cr_ref > 0
        if is_fault || is_desired
            S_ref                   = S;
            S_ref.center_point_type = "Central";
            GEO_ref                 = GEO;
            GEO_ref.nr_fault        = 0;
            GEO_ref.fault_coord_x   = [];
            GEO_ref.fault_coord_y   = [];
            if is_fault
                fprintf('compute_storage_capacity: recomputing Central pressure for reference baseline (Fault mode)...\n');
            else
                fprintf('compute_storage_capacity: recomputing Central pressure for reference baseline (Desired mode)...\n');
            end
            [P_ref, ~] = compute_pressure_buildup(S_ref, R, GEO_ref);
        else
            P_ref = P;
        end

        [OUT.Q_M_each_ref, OUT.V_M_ref, ~, ~, ~, ~, ~, ~] = ...
            capacity_one_mode(P_ref.d_list, P_ref.well_list, P_ref.x_grid_list, P_ref.y_grid_list, P_ref.d_max, ...
                P_ref.p_sup_vec, P_ref.b, P_ref.Q0_vec, ...
                GEO.deltap_cr_ref, R.por, R.thick, S.time_yr, time, R.dens_c, S.maxQ, ...
                'Central', [], []);
        OUT.deltap_cr_ref = GEO.deltap_cr_ref;
        fprintf('compute_storage_capacity: Central det. reference done (Dp_cr_ref = %.2f MPa)\n', ...
            GEO.deltap_cr_ref);
    else
        OUT.Q_M_each_ref  = [];
        OUT.V_M_ref       = [];
        OUT.deltap_cr_ref = NaN;
        warning('compute_storage_capacity: GEO.deltap_cr_ref missing or invalid — reference baseline skipped.');
    end

    % pack output
    OUT.d_list = P.d_list;
    OUT.well_list = P.well_list;
    OUT.d_max = P.d_max;
    OUT.x_grid_list = P.x_grid_list;
    OUT.y_grid_list = P.y_grid_list;
    OUT.distance_matrix = P.distance_matrix;
    OUT.x_length = P.x_length;
    OUT.y_length = P.y_length;
    OUT.well_coords_x = P.well_coords_x;
    OUT.well_coords_y = P.well_coords_y;
    OUT.b = P.b;
    OUT.Q0_vec = P.Q0_vec;
    OUT.p_c = P.p_c;
    OUT.M0 = P.M0;
    OUT.p_sup_vec = P.p_sup_vec;
    OUT.p_sup_2Dgrid = P.p_sup_2Dgrid;
    OUT.Mesh_grid = P.Mesh_grid;
    OUT.p_fault = P.p_fault;
    OUT.fault_coord_x = GEO.fault_coord_x;
    OUT.fault_coord_y = GEO.fault_coord_y;
    OUT.Q_M_each_det = Q_M_each_det;
    OUT.V_M_det = V_M_det;
    OUT.Table_Q_det = Table_Q_det;
    OUT.Table_V_det = Table_V_det;
    OUT.Q_M_each_fault_det = Q_M_each_fault_det;
    OUT.V_M_fault_det = V_M_fault_det;
    OUT.Table_Q_fault_det = Table_Q_fault_det;
    OUT.Table_V_fault_det = Table_V_fault_det;
    OUT.crit_ratio_det = CR.crit_ratio_det;
    OUT.max_ratio_det  = CR.max_ratio_det;
    OUT.p1_fault_det   = CR.p1_fault_det;
    OUT.p2_fault_det   = CR.p2_fault_det;
    OUT.most_critical_fault_table = CR.most_critical_fault_table;
    OUT.max_ratio_P10 = CR.max_ratio_P10;
    OUT.max_ratio_P50 = CR.max_ratio_P50;
    OUT.max_ratio_P90 = CR.max_ratio_P90;
    OUT.Q_M_each_joint_P10 = Q_M_each_joint_P10;  OUT.V_M_joint_P10 = V_M_joint_P10;
    OUT.Q_M_each_joint_P50 = Q_M_each_joint_P50;  OUT.V_M_joint_P50 = V_M_joint_P50;
    OUT.Q_M_each_joint_P90 = Q_M_each_joint_P90;  OUT.V_M_joint_P90 = V_M_joint_P90;
    OUT.Q_M_each_joint_fault_P10 = Q_M_each_joint_fault_P10;  OUT.V_M_joint_fault_P10 = V_M_joint_fault_P10;
    OUT.Q_M_each_joint_fault_P50 = Q_M_each_joint_fault_P50;  OUT.V_M_joint_fault_P50 = V_M_joint_fault_P50;
    OUT.Q_M_each_joint_fault_P90 = Q_M_each_joint_fault_P90;  OUT.V_M_joint_fault_P90 = V_M_joint_fault_P90;
    if ~isempty(p2_fault_mc)
        OUT.p1_fault_mc_P10 = p1_fault_mc(:,:,1);
        OUT.p1_fault_mc_P50 = p1_fault_mc(:,:,2);
        OUT.p1_fault_mc_P90 = p1_fault_mc(:,:,3);
        OUT.p2_fault_mc_P10 = p2_fault_mc(:,:,1);
        OUT.p2_fault_mc_P50 = p2_fault_mc(:,:,2);
        OUT.p2_fault_mc_P90 = p2_fault_mc(:,:,3);
        OUT.mc_idx_fault_P10 = mc_idx_fault(:,:,1);
        OUT.mc_idx_fault_P50 = mc_idx_fault(:,:,2);
        OUT.mc_idx_fault_P90 = mc_idx_fault(:,:,3);
    else
        OUT.p1_fault_mc_P10 = [];  OUT.p1_fault_mc_P50 = [];  OUT.p1_fault_mc_P90 = [];
        OUT.p2_fault_mc_P10 = [];  OUT.p2_fault_mc_P50 = [];  OUT.p2_fault_mc_P90 = [];
        OUT.mc_idx_fault_P10 = []; OUT.mc_idx_fault_P50 = []; OUT.mc_idx_fault_P90 = [];
    end
    OUT.Vmax_fault_samp = Vmax_fault_samp;

    % metadata
    OUT.center_point_type = S.center_point_type;
    OUT.desired_x = S.desired_x;
    OUT.desired_y = S.desired_y;
    OUT.fault_type = S.fault_type;
    OUT.stress_uncertain = S.stress_uncertain;
    OUT.fault_uncertainty = S.fault_uncertainty;
    OUT.n_MC = S.n_MC;
    OUT.n_MC_set = isfield(S,'n_MC_set') && S.n_MC_set;
    OUT.perm_uncertain = MC.perm_active;
    OUT.perm_samples = MC.perm_samples;
    OUT.GEOM = GEOM;
    OUT.no_faults_in_domain = isfield(GEOM,'no_faults_in_domain') && GEOM.no_faults_in_domain;
    if isfield(GEO,'fault_azi');         OUT.fault_azi         = GEO.fault_azi;
    else;                                OUT.fault_azi         = []; end
    if isfield(GEO,'fault_dip');         OUT.fault_dip         = GEO.fault_dip;
    else;                                OUT.fault_dip         = []; end
    if isfield(GEO,'fault_length');      OUT.fault_length      = GEO.fault_length;
    else;                                OUT.fault_length      = []; end
    if isfield(GEO,'deltap_cr_samples'); OUT.deltap_cr_samples = GEO.deltap_cr_samples;
    else;                                OUT.deltap_cr_samples = {}; end
    if isfield(GEO,'deltap_cr_det');     OUT.deltap_cr_det     = GEO.deltap_cr_det;
    else;                                OUT.deltap_cr_det     = []; end

    fprintf('compute_storage_capacity: done [perm_MC=%s | Reference point: %s]\n', ...
        string(logical(MC.perm_active)).replace('true','yes').replace('false','no'), S.center_point_type);
end


function [Q_M_each, V_M, Table_Q, Table_V, ...
          Q_M_each_fault, V_M_fault, Table_Q_fault, Table_V_fault] = ...
    capacity_one_mode(d_list, well_list, x_grid_list, y_grid_list, ...
        d_max, p_sup_vec, b, Q0_vec, ...
        deltap_cr, por, thick, time_yr, time, dens_c, maxQ, ...
        center_point_type, p1_fault, p2_fault)

    nr_dist  = numel(d_list);
    nw       = numel(well_list);
    b_mat    = repmat(b(:),         1, nr_dist);
    q1       = repmat(Q0_vec(:),    1, nr_dist);
    well_mat = repmat(well_list(:), 1, nr_dist);

    p1 = reshape(pad_or_trim(p_sup_vec(:), nw*nr_dist), nw, nr_dist) * 1e6;
    p2 = repmat(deltap_cr * 1e6, size(p1));
    lw_arg = -p2 ./ q1 ./ b_mat .* exp(-p1 ./ q1 ./ b_mat);
    valid  = lw_arg >= -1/exp(1) & lw_arg < 0 & ~isnan(lw_arg);
    q2 = NaN(size(b_mat));
    q2(valid) = real(-p2(valid) ./ b_mat(valid) ./ lambertw(-1, lw_arg(valid)));
    Q_M_each = q2 * 86400*365*dens_c/1e9;
    Q_M_each = apply_constraints(Q_M_each, q2, d_list, well_list, por, thick, time, maxQ);
    V_M      = (Q_M_each .* well_mat) .* time_yr / 1000;
    [Table_Q, Table_V] = make_tables(Q_M_each, V_M, d_list, x_grid_list, y_grid_list);

    [Q_M_each_fault, V_M_fault, Table_Q_fault, Table_V_fault] = ...
        fault_capacity(center_point_type, p1_fault, p2_fault, ...
            q1, b_mat, well_mat, d_list, x_grid_list, y_grid_list, ...
            dens_c, por, thick, time, time_yr, maxQ, nw, nr_dist);
end


function [Q_M_f, V_M_f, TQ_f, TV_f] = ...
    fault_capacity(center_point_type, p1_fault, p2_fault, ...
        q1, b_mat, well_mat, d_list, x_grid_list, y_grid_list, ...
        dens_c, por, thick, time, time_yr, maxQ, nw, nr_dist)

    Q_M_f = [];  V_M_f = [];  TQ_f = table();  TV_f = table();
    if center_point_type ~= "Fault" || isempty(p1_fault) || isempty(p2_fault)
        return;
    end
    if ~isequal(size(p1_fault), size(p2_fault)), return; end

    n   = nw * nr_dist;
    pf1 = reshape(pad_or_trim(p1_fault(:), n), nw, nr_dist);
    pf2 = reshape(pad_or_trim(p2_fault(:), n), nw, nr_dist);
    lw_arg = -pf2 ./ q1 ./ b_mat .* exp(-pf1 ./ q1 ./ b_mat);
    valid  = lw_arg >= -1/exp(1) & lw_arg < 0 & ~isnan(lw_arg);
    q2f = NaN(size(b_mat));
    q2f(valid) = real(-pf2(valid) ./ b_mat(valid) ./ lambertw(-1, lw_arg(valid)));
    Q_M_f = q2f * 86400*365*dens_c/1e9;
    Q_M_f = apply_constraints(Q_M_f, q2f, d_list, x_grid_list .* y_grid_list, por, thick, time, maxQ);
    V_M_f = (Q_M_f .* well_mat) .* time_yr / 1000;
    [TQ_f, TV_f] = make_tables(Q_M_f, V_M_f, d_list, x_grid_list, y_grid_list);
end


function Q_M = apply_constraints(Q_M, q2, d_list, well_list, por, thick, time, maxQ)
    nr_dist = numel(d_list);
    nw      = numel(well_list);
    for dd = 1:nr_dist
        dist = d_list(dd) * 1000;
        if Q_M(1,dd) > 0.9999*maxQ
            Q_M(1,dd) = 0.9999*maxQ;
        end
        for nn = 2:nw
            if isnan(Q_M(nn,dd)), continue; end
            ol = dist^2 * pi * por * thick / 4.001 / time;
            if Q_M(nn,dd) > 0.9999*maxQ || q2(nn,dd) > ol
                Q_M(nn,dd) = min(0.9999*maxQ, ol * 86400*365);
            end
        end
    end
end


function [TQ, TV] = make_tables(Q_M, V_M, d_list, x_grid_list, y_grid_list)
    vn = cellstr(sprintfc('d_%.0f_m', d_list * 1000));
    rn = arrayfun(@(nx,ny) sprintf('%dx%d', nx, ny), ...
        x_grid_list(:), y_grid_list(:), 'UniformOutput', false);
    TQ = array2table(Q_M, 'VariableNames', vn, 'RowNames', rn);
    TQ.Properties.DimensionNames{1} = 'grid_layout';
    TV = array2table(V_M, 'VariableNames', vn, 'RowNames', rn);
    TV.Properties.DimensionNames{1} = 'grid_layout';
end


function v = pad_or_trim(v, n)
    if numel(v) >= n;  v = v(1:n);
    else;              v(end+1:n) = NaN;
    end
end