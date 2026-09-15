function CO2BLOCK2_plots(OUT, R)
% CO2BLOCK_PLOTS  Main plotting function for CO2BLOCK2 results.
%
%   CO2BLOCK2_plots(OUT, R)   called by CO2BLOCK2.m after calculate()
%   CO2BLOCK2_plots()         loads OUT and R from CO2BLOCK2_results.mat
%
%   Window 1  "CO2BLOCK2 Results": tabbed figure with heatmap and contour
%   plots for Q_M and V_M.  Active uncertainty levels determine which tabs
%   appear (Det, Fault-Det, Joint-MC P10/P50/P90).  Stochastic fault mode also
%   adds a fault-length tab.  Duplicate well counts from different grid
%   layouts are collapsed to one row; the best layout per cell is shown on
%   hover.  Excel output (write_tables.m) still includes all grid configs.
%
%   Window 2  "CO2BLOCK2 — Fault Analysis (Interactive)" or
%   "CO2BLOCK2 _ Pressure Map": opened after user picks well scenario and
%   inter-well distance.  Always shows the overpressure map at Q_ref.
%   Additional tabs re-solve ΔP at each scenario rate (Det, Fault-Det,
%   P10/P50/P90), plus a ΔP_cr CDF tab (fault/stress MC) and Perm CDF
%   tab (perm MC).

    if nargin < 1
        loaded = load('CO2BLOCK_results.mat', 'OUT', 'R');
        OUT    = loaded.OUT;
        R      = loaded.R;
        fprintf('Loaded results from CO2BLOCK_results.mat\n');
    end

    %%  global figure defaults 
    set(groot, 'defaultTextInterpreter',          'none');
    set(groot, 'defaultAxesTickLabelInterpreter', 'none');
    set(groot, 'defaultLegendInterpreter',        'none');
    set(groot, 'defaultAxesFontSize',             11);

    % Colour palette
    COL.accent = [0.13 0.47 0.71];   % steel blue
    COL.warm   = [0.84 0.19 0.15];   % brick red
    COL.mid    = [0.20 0.63 0.17];   % forest green
    COL.muted  = [0.45 0.45 0.47];   % mid grey
    COL.bg     = [0.97 0.97 0.975];  % near-white panel background
    COL.amber  = [0.95 0.62 0.12];   % amber

    center_point_type = OUT.center_point_type;
    fault_type        = OUT.fault_type;
    d_list            = OUT.d_list;
    well_list         = OUT.well_list;
    d_max             = OUT.d_max;
    x_grid_list       = OUT.x_grid_list;
    y_grid_list       = OUT.y_grid_list;

    % tick labels show grid layout (e.g. '1x2') so duplicate well counts are clear
    well_labels = arrayfun(@(nx,ny) sprintf('%dx%d', nx, ny), ...
        x_grid_list(:), y_grid_list(:), 'UniformOutput', false);

    %%  convenience flags 
    is_fault    = (center_point_type == "Fault");
    is_desired  = (center_point_type == "Desired");
    any_mc      = isfield(OUT,'Q_M_each_joint_P50') && ~isempty(OUT.Q_M_each_joint_P50);

    %% Window 1- tabbed results figure
    fig_auto = figure('Name', 'CO2BLOCK2 Results', 'Color', 'w', ...
        'Position', [60 60 1300 680]);
    tg_auto = uitabgroup(fig_auto, 'Units','normalized','Position',[0 0 1 1]);

    %%  Tab: Run Summary  (ALWAYS first tab)
    tab = uitab(tg_auto, 'Title', 'Run Summary');
    fig_run_summary(OUT, R, tab, COL);

    %%  Tab: Fault length distribution  (Stochastic fault only) 
    if is_fault && fault_type == "Stochastic"
        fl = [];
        if isfield(OUT,'fault_length') && ~isempty(OUT.fault_length)
            fl = OUT.fault_length(:);
        end
        if ~isempty(fl)
            tab = uitab(tg_auto, 'Title', 'Fault lengths');
            fig_fault_length_dist(fl, COL, tab);
        end
    end

    %%  Tabs: Q_M and V_M - 5 scenarios max 
    stress_unc_active = isfield(OUT,'stress_uncertain') && OUT.stress_uncertain == "yes";

    %% [Det]       Deterministic  - ref perm, mean stress, no fault constraint (always)
    det_QM = safe_out(OUT, 'Q_M_each_ref', safe_out(OUT, 'Q_M_each_det', []));
    det_VM = safe_out(OUT, 'V_M_ref',      safe_out(OUT, 'V_M_det',      []));
    add_QV_tabs(tg_auto, det_QM, det_VM, d_list, well_list, well_labels, d_max, ...
        x_grid_list, y_grid_list, ...
        'Det', ...
        'Deterministic  |  ref perm, mean stress, no fault constraint', COL);

    %% [Fault-Det]  Fault-constrained deterministic (Fault mode only)
    if is_fault
        add_QV_tabs(tg_auto, ...
            safe_out(OUT,'Q_M_each_fault_det',[]), safe_out(OUT,'V_M_fault_det',[]), ...
            d_list, well_list, well_labels, d_max, ...
            x_grid_list, y_grid_list, ...
            'Fault-Det', ...
            'Fault-constrained  |  ref perm, mean stress, most critical fault (deterministic)', COL);
    end

    %% [Det]        Desired-point result (Desired mode only)
    if is_desired
        add_QV_tabs(tg_auto, ...
            safe_out(OUT,'Q_M_each_det',[]), safe_out(OUT,'V_M_det',[]), ...
            d_list, well_list, well_labels, d_max, ...
            x_grid_list, y_grid_list, ...
            'Desired', ...
            'Desired-point result  |  ref perm, mean stress, user (x,y) coordinates', COL);
    end

    %%  [Joint-MC]   P10 / P50 / P90 (any MC active)
    if any_mc
        % choose reservoir or fault-constrained joint depending on mode
        if is_fault
            joint_sets = {
                {'Joint-MC P10', safe_out(OUT,'Q_M_each_joint_fault_P10', safe_out(OUT,'Q_M_each_joint_P10',[])), ...
                              safe_out(OUT,'V_M_joint_fault_P10',      safe_out(OUT,'V_M_joint_P10',[]))}, ...
                {'Joint-MC P50', safe_out(OUT,'Q_M_each_joint_fault_P50', safe_out(OUT,'Q_M_each_joint_P50',[])), ...
                              safe_out(OUT,'V_M_joint_fault_P50',      safe_out(OUT,'V_M_joint_P50',[]))}, ...
                {'Joint-MC P90', safe_out(OUT,'Q_M_each_joint_fault_P90', safe_out(OUT,'Q_M_each_joint_P90',[])), ...
                              safe_out(OUT,'V_M_joint_fault_P90',      safe_out(OUT,'V_M_joint_P90',[]))}};
            joint_subtitle = 'Joint MC  |  per-realization (perm + geomech)  |  fault-constrained';
        else
            joint_sets = {
                {'Joint-MC P10', safe_out(OUT,'Q_M_each_joint_P10',[]), safe_out(OUT,'V_M_joint_P10',[])}, ...
                {'Joint-MC P50', safe_out(OUT,'Q_M_each_joint_P50',[]), safe_out(OUT,'V_M_joint_P50',[])}, ...
                {'Joint-MC P90', safe_out(OUT,'Q_M_each_joint_P90',[]), safe_out(OUT,'V_M_joint_P90',[])}};
            joint_subtitle = 'Joint MC  |  per-realization  |  reservoir-limited';
        end
        for k = 1:3
            lbl = joint_sets{k}{1};
            add_QV_tabs(tg_auto, joint_sets{k}{2}, joint_sets{k}{3}, ...
                d_list, well_list, well_labels, d_max, ...
                x_grid_list, y_grid_list, ...
                lbl, ...
                sprintf('%s  |  %s', joint_subtitle, lbl), COL);
        end
    end

    %% Window 2 - pressure map + fault/perm analysis tabs
    drawnow;
    fprintf('\n=== Pressure Map: Scenario Selection=====\n');
    fprintf('Available well scenarios:\n');
    for i = 1:length(well_list)
        dists_i = OUT.distance_matrix(i, :);
        n_valid = sum(dists_i > 0 & ~isnan(dists_i));
        if n_valid > 0
            fprintf('  [%d]  %s  (%d wells, %d distances)\n', i, well_labels{i}, well_list(i), n_valid);
        else
            fprintf('  [%d]  %s  (%d wells, no valid distances — d_max < dist_min)\n', i, well_labels{i}, well_list(i));
        end
    end

    % let the user retry if they pick a scenario with no valid distances
    while true
        wid = input('Select scenario index: ');
        assert(wid >= 1 && wid <= length(well_list), 'Invalid scenario index.');

        distances = OUT.distance_matrix(wid, :);
        valid_idx = find(distances > 0 & ~isnan(distances));
        if ~isempty(valid_idx)
            break;
        end
        fprintf('  Scenario %s has no valid distances (d_max < dist_min). Pick another.\n', well_labels{wid});
    end

    fprintf('Available inter-well distances:\n');
    for i = valid_idx
        fprintf('  [%d]  %.1f km\n', i, distances(i)/1000);
    end
    did = input('Select distance index: ');
    assert(ismember(did, valid_idx), 'Invalid distance index.');

    % Open Window 2 with a title that reflects the active mode
    if is_fault
        win2_name = 'CO2BLOCK2 — Fault Analysis (Interactive)';
    else
        win2_name = 'CO2BLOCK2 — Pressure Map';
    end
    fig_int = figure('Name', win2_name, 'Color', 'w', 'Position', [100 50 1200 720]);
    tg_int  = uitabgroup(fig_int, 'Units','normalized','Position',[0 0 1 1]);

    %  Tab: Overpressure map at Q_ref  (ALWAYS — all modes) ──────────
    fprintf('Computing pressure map at Q_ref ...\n');
    [p_map_ref, X, Y] = compute_pressure_map(OUT.GEOM, R, wid, did);
    tab = uitab(tg_int, 'Title', 'Overpressure at Q_{ref}');
    fig_overpressure_fault_map(p_map_ref, X, Y, OUT, wid, did, COL, tab, struct('show_crit', false));

    % ── Tabs: ΔP re-solved at each Q_M scenario rate ─────────────────
    % Pressure is NOT proportional to Q — Q enters the plume radius
    % (csi/psi) as well as the amplitude — so each map is solved from
    % scratch at its own rate rather than rescaled from p_map_ref.
    add_scaled_pressure_tabs(tg_int, p_map_ref, X, Y, OUT, R, wid, did, COL, is_fault, any_mc);

    if isfield(OUT,'Q_M_each_fault_det') && ~isempty(OUT.Q_M_each_fault_det) && ...
       isfield(OUT,'Q_M_each_det')        && ~isempty(OUT.Q_M_each_det)
        Q_max_pw = OUT.Q_M_each_fault_det(wid, did);
        Q_ref_pw = OUT.Q_M_each_det(wid, did);
        if isfinite(Q_max_pw) && isfinite(Q_ref_pw) && Q_ref_pw > 0
            fprintf('  Q_max/Q_ref scale = %.4f\n', Q_max_pw / Q_ref_pw);
        end
    end

    %  Tab: ΔP_cr CDF  (Fault mode OR non-fault stress uncertainty) ──
    if is_fault || stress_unc_active
        tab = uitab(tg_int, 'Title', 'ΔP_cr CDF');
        fig_deltap_cr_cdf(OUT, COL, wid, did, tab);
    end

    %  Tab: Perm CDF  (any mode, whenever perm MC is active) =
    if any_mc && isfield(OUT,'perm_samples') && ~isempty(OUT.perm_samples)
        tab = uitab(tg_int, 'Title', 'Perm CDF');
        fig_perm_cdf(OUT, COL, tab);
    end

    fprintf('\nAll CO2BLOCK2 plots generated.\n');
end


%% safe_out
function v = safe_out(s, fname, default)
    if isfield(s, fname) && ~isempty(s.(fname))
        v = s.(fname);
    else
        v = default;
    end
end


%% aggregate_by_wellcount  - collapses duplicate well counts into one row
function [QM_agg, VM_agg, wl_agg, labels_agg, best_layout] = aggregate_by_wellcount( ...
        QM, VM, well_list, x_grid_list, y_grid_list)

    [unique_wc, ~, ic] = unique(well_list(:));
    nw_agg = numel(unique_wc);
    nd     = size(QM, 2);

    QM_agg      = NaN(nw_agg, nd);
    VM_agg      = NaN(nw_agg, nd);
    wl_agg      = unique_wc(:);
    labels_agg  = cell(nw_agg, 1);
    best_layout = cell(nw_agg, nd);

    for g = 1:nw_agg
        idx = find(ic == g);
        labels_agg{g} = sprintf('%d', unique_wc(g));

        if numel(idx) == 1
            QM_agg(g, :) = QM(idx, :);
            if ~isempty(VM), VM_agg(g, :) = VM(idx, :); end
            layout_str = sprintf('%dx%d', x_grid_list(idx), y_grid_list(idx));
            for di = 1:nd
                best_layout{g, di} = layout_str;
            end
        else
            for di = 1:nd
                col_vals = QM(idx, di);
                finite_mask_col = isfinite(col_vals);
                if any(finite_mask_col)
                    [~, best_k] = max(col_vals);
                else
                    best_k = 1;
                end
                best_row = idx(best_k);
                QM_agg(g, di) = QM(best_row, di);
                if ~isempty(VM), VM_agg(g, di) = VM(best_row, di); end
                best_layout{g, di} = sprintf('%dx%d', ...
                    x_grid_list(best_row), y_grid_list(best_row));
            end
        end
    end

    if isempty(VM)
        VM_agg = [];
    end
end


%% add_QV_tabs
function add_QV_tabs(tg, QM, VM, d_list, well_list, well_labels, d_max, ...
                     x_grid_list, y_grid_list, tab_lbl, chart_subtitle, COL)
    if isempty(QM), return; end

    [QM_agg, VM_agg, wl_agg, labels_agg, best_layout] = ...
        aggregate_by_wellcount(QM, VM, well_list, x_grid_list, y_grid_list);

    tab = uitab(tg, 'Title', sprintf('Q_M [%s] heatmap', tab_lbl));
    draw_heatmap_tab(tab, real(QM_agg), d_list, wl_agg, labels_agg, d_max, chart_subtitle, COL, 'Q_M  (Mt/yr)', best_layout);

    tab = uitab(tg, 'Title', sprintf('Q_M [%s] contour', tab_lbl));
    draw_contour_tab(tab, real(QM_agg), d_list, wl_agg, labels_agg, d_max, chart_subtitle, ...
        'Q_M  (Mt/yr)', false, COL);

    if ~isempty(VM_agg) && any(isfinite(VM_agg(:)))
        tab = uitab(tg, 'Title', sprintf('V_M [%s] heatmap', tab_lbl));
        draw_heatmap_tab(tab, real(VM_agg), d_list, wl_agg, labels_agg, d_max, chart_subtitle, COL, 'V_M  (Gt)', best_layout);

        tab = uitab(tg, 'Title', sprintf('V_M [%s] contour', tab_lbl));
        draw_contour_tab(tab, real(VM_agg), d_list, wl_agg, labels_agg, d_max, chart_subtitle, ...
            'V_M  (Gt)', true, COL);
    end
end


%% fig_fault_length_dist
function fig_fault_length_dist(fault_length, COL, parent)

    L_km = fault_length / 1000;
    n    = numel(L_km);

    ax = axes('Parent', parent, 'Position', [0.12 0.13 0.82 0.74]); %#ok<LAXES>
    hold(ax, 'on');

    histogram(ax, L_km, 'NumBins', 30, 'FaceColor', COL.accent, ...
        'EdgeColor', 'w', 'LineWidth', 0.4, 'FaceAlpha', 0.85);

    xline(mean(L_km),   '--', 'Color', COL.warm, 'LineWidth', 2.0, ...
        'Label', sprintf('Mean = %.2f km', mean(L_km)), ...
        'LabelOrientation', 'horizontal', 'FontSize', 9, ...
        'LabelVerticalAlignment', 'top');
    xline(median(L_km), ':',  'Color', COL.mid,  'LineWidth', 2.0, ...
        'Label', sprintf('Median = %.2f km', median(L_km)), ...
        'LabelOrientation', 'horizontal', 'FontSize', 9, ...
        'LabelVerticalAlignment', 'bottom');

    xlabel(ax, 'Fault Length  (km)', 'FontSize', 12);
    ylabel(ax, 'Number of Faults',   'FontSize', 12);
    title(ax, ...
        sprintf('Stochastic Fault Population — Length Distribution  (n = %d faults)', n), ...
        'FontSize', 13, 'FontWeight', 'bold');
    set(ax, 'Color', COL.bg, 'Box', 'on'); grid on; grid minor;

    pcts = [10 50 90];
    pclr = {COL.mid, COL.amber, COL.warm};
    for ii = 1:3
        pv = prctile(L_km, pcts(ii));
        xline(pv, '-', 'Color', [pclr{ii} 0.5], 'LineWidth', 1.0, ...
            'Label', sprintf('P%d=%.2f', pcts(ii), pv), ...
            'FontSize', 8, 'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom', 'HandleVisibility','off');
    end

    stats_str = sprintf('n = %d\nMin: %.2f km\nMax: %.2f km\nMean: %.2f km\nStd: %.2f km', ...
        n, min(L_km), max(L_km), mean(L_km), std(L_km));
    text(ax, 0.97, 0.97, stats_str, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'BackgroundColor', 'w', ...
        'EdgeColor', [0.70 0.70 0.70], 'Margin', 4);
    hold(ax, 'off');
end


%% fig_deltap_cr_cdf
function fig_deltap_cr_cdf(OUT, COL, wid, did, parent)

    ax = axes('Parent', parent, 'Position', [0.12 0.13 0.82 0.74]); %#ok<LAXES>
    hold(ax, 'on');

    have_full = isfield(OUT,'deltap_cr_samples') && ...
                ~isempty(OUT.deltap_cr_samples) && ...
                numel(OUT.deltap_cr_samples) >= 1;

    if have_full
        if iscell(OUT.deltap_cr_samples)
            try
                raw_r = OUT.crit_ratio_det{wid}{did}(:);
                [~, idx_c] = max(raw_r);
            catch
                idx_c = 1;
            end
            idx_c      = min(idx_c, numel(OUT.deltap_cr_samples));
            samples_mc = OUT.deltap_cr_samples{idx_c};
        else
            samples_mc = OUT.deltap_cr_samples(:);
        end
        samples_mc = samples_mc(isfinite(samples_mc) & samples_mc > 0);

        vals   = sort(samples_mc);
        n      = numel(vals);
        cum_p  = (1:n) / n * 100;

        plot(ax, vals, cum_p, '-', 'Color', COL.accent, 'LineWidth', 2.5, ...
            'DisplayName', sprintf('\\DeltaP_{cr} — most critical fault  (n = %d MC realisations)', n));

        pcts     = [10 50 90];
        pct_clrs = {COL.mid, COL.amber, COL.warm};
        pct_lbls = {'P10 — conservative estimate (90% of realisations give higher Δp_cr)', ...
                    'P50 — median estimate', ...
                    'P90 — optimistic estimate (only 10% of realisations give higher Δp_cr)'};
        for ii = 1:3
            pv = prctile(samples_mc, pcts(ii));
            xline(ax, pv, '--', 'Color', pct_clrs{ii}, 'LineWidth', 1.6, ...
                'HandleVisibility','off');
            scatter(ax, pv, pcts(ii), 80, pct_clrs{ii}, 'filled', ...
                'MarkerEdgeColor','w', 'LineWidth',0.5, ...
                'DisplayName', sprintf('P%d = %.2f MPa  — %s', pcts(ii), pv, pct_lbls{ii}));
            text(ax, pv, pcts(ii)+2, sprintf('  P%d=%.2f MPa', pcts(ii), pv), ...
                'FontSize', 8, 'Color', pct_clrs{ii});
        end
        if iscell(OUT.deltap_cr_samples)
            title_note = sprintf('Most critical fault in scenario: %d wells, %.1f km spacing', ...
                OUT.well_list(wid), OUT.distance_matrix(wid,did)/1000);
        else
            title_note = sprintf('Central/Desired ref point — stress uncertainty  (n = %d MC realisations)', ...
                numel(samples_mc));
        end

    else
        fprintf(['Fig 2c: deltap_cr_samples not in OUT — using p2_fault across all\n' ...
            'scenarios as proxy.  See function header for how to enable full MC CDF.\n']);

        if isfield(OUT,'p2_fault_det') && ~isempty(OUT.p2_fault_det)
            sets  = {OUT.p2_fault_det};
            names = {'\DeltaP_{cr} — deterministic (mean stress)'};
            clrs  = {COL.accent};
        else
            sets = {}; names = {}; clrs = {};
        end

        for k = 1:numel(sets)
            if isempty(sets{k}), continue; end
            vals  = sets{k}(:) / 1e6;
            vals  = sort(vals(isfinite(vals) & vals > 0));
            if isempty(vals), continue; end
            n     = numel(vals);
            cum_p = (1:n) / n * 100;
            plot(ax, vals, cum_p, '-', 'Color', clrs{k}, 'LineWidth', 2.2, ...
                'DisplayName', names{k});
            for pp = [10 50 90]
                pv = prctile(vals, pp);
                scatter(ax, pv, pp, 55, clrs{k}, 'filled', ...
                    'MarkerEdgeColor','w','HandleVisibility','off');
                text(ax, pv, pp+2, sprintf('  P%d=%.1f', pp, pv), ...
                    'FontSize', 7, 'Color', clrs{k});
            end
        end
        title_note = 'Distribution across all (scenario, distance) combinations';
    end

    xlabel(ax, 'Critical Reactivation Overpressure  ΔP_cr  (MPa)', 'FontSize', 12);
    ylabel(ax, 'Cumulative Frequency  (%)', 'FontSize', 12);
    title(ax, 'CDF of Critical Overpressure  ΔP_cr', ...
        'FontSize', 13, 'FontWeight', 'bold');
    subtitle(ax, [title_note '  |  ' ...
        'ΔP_cr = excess reservoir pressure that would just reactivate the fault'], ...
        'FontSize', 9);
    legend(ax, 'Location', 'southeast', 'FontSize', 9, 'Interpreter', 'tex');
    ylim(ax, [0 100]); grid(ax,'on'); grid(ax,'minor');
    set(ax, 'Color', COL.bg, 'Box', 'on');
    hold(ax, 'off');
end


%% fig_perm_cdf
function fig_perm_cdf(OUT, COL, parent)

    if ~isfield(OUT,'perm_samples') || isempty(OUT.perm_samples)
        ax = axes('Parent', parent); %#ok<LAXES>
        text(ax, 0.5, 0.6, ...
            {'Permeability CDF requires OUT.perm\_samples to be saved.', ...
             '', ...
             'Add this line to compute\_storage\_capacity.m', ...
             'just before the final fprintf(...):  ', ...
             '', ...
             '   OUT.perm\_samples = MC.perm\_samples;', ...
             '', ...
             'Then re-run CO2BLOCK2.m to generate this figure.'}, ...
            'Units','normalized','HorizontalAlignment','center', ...
            'VerticalAlignment','middle','FontSize',11,'FontName','Monospaced');
        title(ax,'Permeability CDF — data not yet saved','FontSize',12);
        axis(ax,'off');
        return;
    end

    perm_mD = OUT.perm_samples(:) / 1e-15;
    perm_mD = sort(perm_mD(perm_mD > 0));
    n       = numel(perm_mD);
    cum_p   = (1:n) / n * 100;

    ax = axes('Parent', parent, 'Position', [0.13 0.13 0.81 0.74]); %#ok<LAXES>
    hold(ax,'on');

    plot(ax, perm_mD, cum_p, '-', 'Color', COL.accent, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Permeability CDF  (n = %d realisations)', n));

    pcts     = [10 50 90];
    pct_clrs = {COL.mid, COL.amber, COL.warm};
    pct_lbls = {'P10 — low permeability (conservative capacity estimate)', ...
                'P50 — median permeability', ...
                'P90 — high permeability (optimistic capacity estimate)'};
    for ii = 1:3
        pv = prctile(perm_mD, pcts(ii));
        xline(ax, pv, '--', 'Color', pct_clrs{ii}, 'LineWidth', 1.6, ...
            'HandleVisibility','off');
        scatter(ax, pv, pcts(ii), 80, pct_clrs{ii}, 'filled', ...
            'MarkerEdgeColor','w','LineWidth',0.5, ...
            'DisplayName', sprintf('P%d = %.3g mD  — %s', pcts(ii), pv, pct_lbls{ii}));
        text(ax, pv*1.02, pcts(ii)+2, sprintf('P%d=%.3g mD', pcts(ii), pv), ...
            'FontSize', 8, 'Color', pct_clrs{ii});
    end

    if isfield(OUT,'M0') && ~isempty(OUT.M0)
        ref_mD = OUT.M0 * 100;
        if isfinite(ref_mD) && ref_mD > 0
            xline(ax, ref_mD, '-k', 'LineWidth', 1.5, ...
                'Label', sprintf('Reference perm = %.3g mD', ref_mD), ...
                'FontSize', 9, 'LabelVerticalAlignment', 'bottom');
        end
    end

    set(ax, 'XScale', 'log', 'Color', COL.bg, 'Box', 'on');
    xlabel(ax, 'Permeability  (mD, log scale)', 'FontSize', 12);
    ylabel(ax, 'Realisations with permeability below this value  (%)', 'FontSize', 12);
    title(ax, 'CDF of Permeability — Monte Carlo Realisations', ...
        'FontSize', 13, 'FontWeight', 'bold');
    subtitle(ax, ['Each realisation draws a permeability from the log-normal distribution ' ...
        'and computes Q_M and V_M independently.  ' ...
        'P10/P50/P90 define the low/median/high permeability scenarios and their capacity estimates.'], ...
        'FontSize', 9);
    legend(ax, 'Location', 'southeast', 'FontSize', 9);
    ylim(ax,[0 100]); grid(ax,'on'); grid(ax,'minor');
    hold(ax,'off');
end


%% fig_overpressure_fault_map
% opts fields (all optional): title_label, dc_ref, dc_vec, p_fault_scen, dc_note.
% dc_vec and p_fault_scen must come from the same realization.
% Domain centre (x_length/2, y_length/2) always marked as a star.
% Fault segment data tips also show Criticality R on hover.
function fig_overpressure_fault_map(p_map, X, Y, OUT, wid, did, COL, parent, opts)

    %%  optional display overrides 
    if nargin < 9 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'title_label'),  opts.title_label  = 'Q_{ref}'; end
    if ~isfield(opts,'dc_ref'),       opts.dc_ref       = [];        end
    if ~isfield(opts,'dc_vec'),       opts.dc_vec       = [];        end
    if ~isfield(opts,'p_fault_scen'), opts.p_fault_scen = [];        end
    if ~isfield(opts,'dc_note'),      opts.dc_note      = '';        end
        if ~isfield(opts,'show_crit'),    opts.show_crit    = true;      end

    is_fault_mode = (OUT.center_point_type == "Fault") && ...
                    isfield(OUT.GEOM,'nr_fault') && (OUT.GEOM.nr_fault > 0);

    if is_fault_mode && opts.show_crit
        ax_pos = [0.08 0.12 0.58 0.78];
    else
        ax_pos = [0.08 0.12 0.70 0.78];
    end

    ax_main = axes('Parent', parent, 'Position', ax_pos); %#ok<LAXES>
    hold(ax_main, 'on');

    %%  background overpressure contour 
    [~, hcf] = contourf(ax_main, X/1000, Y/1000, p_map, 25, 'LineColor', 'none');
    hcf.HandleVisibility = 'off';
    colormap(ax_main, make_pressure_cmap(256));

    % ── colorbar label (single line — MATLAB rotates it, newlines render poorly) ──
    cb1_str = 'Overpressure   ΔP   (MPa)';

    %%  injection wells (always shown)
    wx = OUT.well_coords_x{wid,did}(:) / 1000;
    wy = OUT.well_coords_y{wid,did}(:) / 1000;
    scatter(ax_main, wx, wy, 90, COL.muted, '^', 'filled', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.6, ...
        'DisplayName', 'Injection wells');

    %%  Reference point marker - ONE marker matching center_point_type 
    % Fault mode  : most critical fault = reference point; marked below.
    % Central mode: geometric centre of reservoir domain (x_length/2, y_length/2).
    % Desired mode: user-specified (x,y) from General_settings.xlsx.
    if OUT.center_point_type == "Central"
        cx_km = max(X(:)) / 2 / 1000;
        cy_km = max(Y(:)) / 2 / 1000;
        hs = scatter(ax_main, cx_km, cy_km, ...
            220, COL.accent, 'h', 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Reference point  (domain centre: %.1f, %.1f km)', cx_km, cy_km));
        hs.DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('X  (km)', cx_km), ...
            dataTipTextRow('Y  (km)', cy_km)];

    elseif OUT.center_point_type == "Desired"
        xd_km = []; yd_km = [];
        if isfield(OUT,'desired_x') && isfield(OUT,'desired_y')
            xd_km = OUT.desired_x;  yd_km = OUT.desired_y;
        elseif isfield(OUT,'GEOM') && isfield(OUT.GEOM,'desired_x')
            xd_km = OUT.GEOM.desired_x;  yd_km = OUT.GEOM.desired_y;
        end
        if ~isempty(xd_km) && isfinite(xd_km)
            hs = scatter(ax_main, xd_km, yd_km, ...
                260, COL.amber, 'p', 'filled', ...
                'MarkerEdgeColor', 'w', 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Reference point  (desired: %.1f, %.1f km)', xd_km, yd_km));
            hs.DataTipTemplate.DataTipRows = [ ...
                dataTipTextRow('X  (km)', xd_km), ...
                dataTipTextRow('Y  (km)', yd_km)];
        end
    end

    if is_fault_mode
        %%  fault geometry 
        nr_fault = OUT.GEOM.nr_fault;
        fx = OUT.fault_coord_x / 1000;
        fy = OUT.fault_coord_y / 1000;

        fazi  = zeros(nr_fault, 1);
        fl_km = 0.5 * ones(nr_fault, 1);
        if isfield(OUT,'fault_azi')    && ~isempty(OUT.fault_azi)
            fazi  = OUT.fault_azi(:);
        end
        if isfield(OUT,'fault_length') && ~isempty(OUT.fault_length)
            fl_km = OUT.fault_length(:) / 1000;
        end

        if opts.show_crit
            %%  criticality ratios 
            crit_arr = get_crit_ratio(OUT, wid, did);   % default: Q_ref, deterministic

            % R = ΔP at this rate / Δp_cr of the SAME realization.  Scaling the
            % numerator to the scenario rate while leaving the denominator
            % deterministic used to let R exceed 1 spuriously.
            if ~isempty(opts.p_fault_scen) && ~isempty(opts.dc_vec)
                n   = min(numel(opts.p_fault_scen), numel(opts.dc_vec));
                dcv = opts.dc_vec(1:n);
                dcv(dcv <= 0) = NaN;
                crit_arr = NaN(nr_fault, 1);
                crit_arr(1:n) = opts.p_fault_scen(1:n) ./ dcv;
            end

            cmap_fault = make_fault_cmap(256);

            if ~isempty(crit_arr) && any(isfinite(crit_arr))
                cmin = max(0, min(crit_arr(isfinite(crit_arr))));
                cmax = max(cmin + 0.01, max(crit_arr(isfinite(crit_arr))));
            else
                cmin = 0; cmax = 1;
            end
            if isempty(crit_arr), crit_arr = NaN(nr_fault, 1); end

            %% draw fault segments coloured by R 
            % Each line handle is captured so a DataTipTemplate can show
            % Criticality R alongside X and Y on hover/click.
            for f = 1:nr_fault
                dx = (fl_km(f)/2) * sind(fazi(f));
                dy = (fl_km(f)/2) * cosd(fazi(f));
                xs = [fx(f)-dx, fx(f)+dx];
                ys = [fy(f)-dy, fy(f)+dy];
                if isfinite(crit_arr(f))
                    t      = min(1, max(0, (crit_arr(f)-cmin) / max(cmax-cmin, 1e-9)));
                    ci     = max(1, round(t*255)+1);
                    fc_rgb = cmap_fault(ci,:);
                    lw     = 1.4 + 2.8*t;
                else
                    fc_rgb = [0.65 0.65 0.65];
                    lw     = 1.2;
                end
                hl = plot(ax_main, xs, ys, '-', 'Color', fc_rgb, 'LineWidth', lw, ...
                    'HandleVisibility', 'off');
                if isfinite(crit_arr(f))
                    hl.DataTipTemplate.DataTipRows(end+1) = ...
                        dataTipTextRow('Criticality R', repmat(crit_arr(f), 1, 2));
                end
            end

            %%  legend line for fault (single neutral gray line) 
            plot(ax_main, NaN, NaN, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 2.5, ...
                'DisplayName', 'Fault');

            %%  most critical fault marker 
            if any(isfinite(crit_arr))
                [max_cr, idx_c] = max(crit_arr);
                hs_f = scatter(ax_main, fx(idx_c), fy(idx_c), 220, COL.warm, 'p', 'filled', ...
                    'LineWidth', 0.9, ...
                    'DisplayName', sprintf('Reference point  |  Most critical fault  (R = %.2f)', max_cr));
                hs_f.DataTipTemplate.DataTipRows = [ ...
                    dataTipTextRow('X  (km)',           fx(idx_c)), ...
                    dataTipTextRow('Y  (km)',           fy(idx_c)), ...
                    dataTipTextRow('Criticality R',     max_cr)];
            end

        else
            %%  faults shown in gray (no criticality coloring) 
            for f = 1:nr_fault
                dx = (fl_km(f)/2) * sind(fazi(f));
                dy = (fl_km(f)/2) * cosd(fazi(f));
                xs = [fx(f)-dx, fx(f)+dx];
                ys = [fy(f)-dy, fy(f)+dy];
                plot(ax_main, xs, ys, '-', 'Color', [0.50 0.50 0.50], 'LineWidth', 1.4, ...
                    'HandleVisibility', 'off');
            end
            plot(ax_main, NaN, NaN, '-', 'Color', [0.50 0.50 0.50], 'LineWidth', 2.5, ...
                'DisplayName', 'Fault');
        end  % show_crit

        %  title uses scenario label ─────────────────────────────────
        title_str = sprintf('ΔP at %s  |  Fault Criticality  |  %d Wells  |  %.1f km', ...
            opts.title_label, OUT.well_list(wid), OUT.distance_matrix(wid,did)/1000);

        unc_parts = {};
        if isfield(OUT,'stress_uncertain') && OUT.stress_uncertain == "yes"
            unc_parts{end+1} = 'stress MC';
        end
        if isfield(OUT,'fault_uncertainty') && OUT.fault_uncertainty == "yes"
            unc_parts{end+1} = 'fault geometry MC';
        end
        if isfield(OUT,'perm_uncertain') && OUT.perm_uncertain
            unc_parts{end+1} = 'perm MC';
        end
        if isempty(unc_parts)
            unc_str = 'Deterministic (no uncertainty)';
        else
            unc_str = ['Uncertainty: ' strjoin(unc_parts, ' + ')];
        end

        % only mention Q_ref on the base tab; drop it when a specific Q is shown
        if strcmp(opts.title_label, 'Q_{ref}')
            sub_str = [unc_str, '  |  Q_{ref} = reservoir-limit injection rate  |  ref perm, mean stress'];
            if isfield(OUT,'max_ratio_P10') && ~isempty(OUT.max_ratio_P10)
                sub_str = [sub_str, sprintf(['  |  R_{max} at Q_{ref}: ' ...
                    'P10 = %.2f, P50 = %.2f, P90 = %.2f'], ...
                    OUT.max_ratio_P10(wid,did), OUT.max_ratio_P50(wid,did), ...
                    OUT.max_ratio_P90(wid,did))];
            end
        else
            sub_str = [unc_str, '  |  ref perm'];
            if ~isempty(opts.dc_note)
                sub_str = [sub_str, '  |  ', opts.dc_note];
            end
            if isfield(OUT,'perm_uncertain') && OUT.perm_uncertain
                sub_str = [sub_str, '  |  NOTE: perm MC active — map drawn at reference perm'];
            end
        end

        if opts.show_crit
            %% secondary colorbar: criticality ratio
            ax2 = axes('Parent', parent, 'Position', [0.01 0.01 0.01 0.01], ...
                'Visible', 'off'); %#ok<LAXES>
            colormap(ax2, cmap_fault);
            cb2 = colorbar(ax2);
            cb2.Ticks      = linspace(0, 1, 6);
            cb2.TickLabels = arrayfun(@(t) sprintf('%.2f', cmin + t*(cmax-cmin)), ...
                cb2.Ticks, 'UniformOutput', false);
            cb2.Label.String   = 'Criticality Ratio   R = ΔP_fault / ΔP_cr';
            cb2.Label.FontSize = 10;
            caxis(ax2, [0 1]);
            rt = min(1, max(0, (1-cmin) / max(cmax-cmin, 1e-9)));
            if rt > 0.02 && rt < 0.98
                hold(ax2,'on');
                plot(ax2, [0 1], [rt rt], '--w', 'LineWidth', 1.2);
                hold(ax2,'off');
            end

        end  % show_crit cb2
    else
        %%  Reservoir / Desired mode 
        title_str = sprintf('Overpressure at %s  |  %d Wells  |  %.1f km spacing', ...
            opts.title_label, OUT.well_list(wid), OUT.distance_matrix(wid,did)/1000);
        unc_parts = {};
        if isfield(OUT,'stress_uncertain') && OUT.stress_uncertain == "yes"
            unc_parts{end+1} = 'stress MC';
        end
        if isfield(OUT,'perm_uncertain') && OUT.perm_uncertain
            unc_parts{end+1} = 'perm MC';
        end
        if isempty(unc_parts)
            sub_str = sprintf('Reference point: %s  |  Deterministic (no uncertainty)', ...
                OUT.center_point_type);
        else
            sub_str = sprintf('Reference point: %s  |  Uncertainty: %s', ...
                OUT.center_point_type, strjoin(unc_parts,' + '));
        end
    end

    daspect(ax_main, [1 1 1]);
    xlim(ax_main, [0, max(X(:))/1000]);
    ylim(ax_main, [0, max(Y(:))/1000]);

    cb_gap  = 0.015;
    cb_w    = 0.022;
    cb2_gap = 0.055;   % extra room for the rotated pressure label
    cb1_x   = ax_pos(1) + ax_pos(3) + cb_gap;
    cb1 = colorbar(ax_main, 'Position', [cb1_x, ax_pos(2), cb_w, ax_pos(4)]);
    cb1.Label.String   = cb1_str;
    cb1.Label.FontSize = 10;
    if exist('cb2','var')
        cb2.Position = [cb1_x + cb_w + cb2_gap, ax_pos(2), cb_w, ax_pos(4)];
    end

    xlabel(ax_main, 'Easting  (km)',  'FontSize', 12);
    ylabel(ax_main, 'Northing  (km)', 'FontSize', 12);
    title(ax_main,   title_str, 'FontSize', 12, 'FontWeight', 'bold');
    subtitle(ax_main, sub_str,  'FontSize',  9);
    legend(ax_main, 'Location', 'best', 'FontSize', 9);
    set(ax_main, 'Box', 'on', 'Color', COL.bg);
    hold(ax_main, 'off');
end


%% add_scaled_pressure_tabs
% Solves the pressure field at each scenario rate and opens a tab in Window 2.
% Each tab pairs ΔP at its own rate with Δp_cr from the same source (det or
% matched MC realization) so R stays consistent. On the Det tab R > 1 is
% expected — that rate is the reservoir limit, fault constraint not applied.
% For Joint-MC P10/P50/P90, fault pressures are scaled by k_ref/k_star to
% account for the realization perm; the background map stays at ref perm.
% Scenarios: Det, Fault-Det, Joint-MC P10/P50/P90 (when data exists and Q > 0).
function add_scaled_pressure_tabs(tg, p_map_ref, X, Y, OUT, R, wid, did, COL, is_fault, any_mc) %#ok<INUSL>

    % {tab label, Q field, Δp_cr field for the reference fault, percentile slot}
    pairs = {};
    pairs{end+1} = {'Det',       'Q_M_each_ref',       'p2_fault_det', 0};
    if is_fault
        pairs{end+1} = {'Fault-Det', 'Q_M_each_fault_det', 'p2_fault_det', 0};
        if any_mc
            pairs{end+1} = {'Joint-MC P10', 'Q_M_each_joint_fault_P10', 'p2_fault_mc_P10', 1};
            pairs{end+1} = {'Joint-MC P50', 'Q_M_each_joint_fault_P50', 'p2_fault_mc_P50', 2};
            pairs{end+1} = {'Joint-MC P90', 'Q_M_each_joint_fault_P90', 'p2_fault_mc_P90', 3};
        end
    elseif any_mc
        pairs{end+1} = {'Joint-MC P10', 'Q_M_each_joint_P10', '', 0};
        pairs{end+1} = {'Joint-MC P50', 'Q_M_each_joint_P50', '', 0};
        pairs{end+1} = {'Joint-MC P90', 'Q_M_each_joint_P90', '', 0};
    end

    for k = 1:numel(pairs)
        lbl = pairs{k}{1};
        fld = pairs{k}{2};
        dcf = pairs{k}{3};
        ip  = pairs{k}{4};
        if ~isfield(OUT, fld) || isempty(OUT.(fld)), continue; end
        qv = OUT.(fld)(wid, did);              % Mt/yr per well
        if ~isfinite(qv) || qv <= 0, continue; end

        % Re-solve the pressure field at this rate.  Q also enters the plume
        % radius (csi/psi), so p is NOT proportional to Q.
        p_scen = compute_pressure_map(OUT.GEOM, R, wid, did, qv);

        sc_opts = struct();
        if ip == 0
            sc_opts.dc_note = 'ΔP_{cr} deterministic (mean stress, nominal fault geometry)';
        else
            sc_opts.dc_note = 'ΔP_{cr} and R from the MC realization matching this rate  |  R corrected to realization perm';
        end
        if strcmp(lbl, 'Det')
            sc_opts.title_label = sprintf(['Det — reservoir limit, fault constraint NOT ' ...
                'applied  (%.3g Mt/yr/well)'], qv);
            sc_opts.dc_note = [sc_opts.dc_note '  |  R > 1 here is expected'];
            sc_opts.show_crit = false;   % faults in gray; no criticality colorbar
        else
            sc_opts.title_label = sprintf('%s  (%.3g Mt/yr/well)', lbl, qv);
        end

        % Δp_cr from the SAME realization that produced this rate
        if ~isempty(dcf) && isfield(OUT, dcf) && ~isempty(OUT.(dcf))
            sc_opts.dc_ref = OUT.(dcf)(wid, did) / 1e6;      % MPa
        end
        sc_opts.dc_vec = get_dc_vector(OUT, wid, did, ip);   % MPa, per fault

        if isfield(OUT,'fault_coord_x') && ~isempty(OUT.fault_coord_x)
            sc_opts.p_fault_scen = pressure_at_points(OUT.GEOM, R, wid, did, qv, ...
                OUT.fault_coord_x(:), OUT.fault_coord_y(:));
        end

        % pressure_at_points always uses reference perm.  For Joint-MC P10/P50/P90
        % tabs the matched realization may have perm != ref perm, making
        % p_fault_scen too large or too small by factor k_ref/k_star.
        % Scale only the fault pressures used for R coloring — the background
        % pressure colormap (p_scen) stays at ref perm for visual comparability.
        if ip > 0 && isfield(OUT,'perm_samples') && ~isempty(OUT.perm_samples)
            fld_idx = {'mc_idx_fault_P10','mc_idx_fault_P50','mc_idx_fault_P90'};
            fld_ip  = fld_idx{ip};
            if isfield(OUT, fld_ip) && ~isempty(OUT.(fld_ip))
                i_star = OUT.(fld_ip)(wid, did);
                if isfinite(i_star) && i_star >= 1 && i_star <= numel(OUT.perm_samples)
                    k_star = OUT.perm_samples(i_star);
                    if k_star > 0 && isfinite(k_star) && isfield(sc_opts,'p_fault_scen')
                        sc_opts.p_fault_scen = sc_opts.p_fault_scen * (R.perm / k_star);
                    end
                end
            end
        end

        tab = uitab(tg, 'Title', sprintf('ΔP [%s]', lbl));
        fig_overpressure_fault_map(p_scen, X, Y, OUT, wid, did, COL, tab, sc_opts);
    end
end


%% get_dc_vector
function dc_vec = get_dc_vector(OUT, wid, did, ip)
% Per-fault Δp_cr [MPa] for the realization behind this tab.
%   ip = 0        -> deterministic (mean stress, nominal geometry)
%   ip = 1/2/3    -> the MC draw matched to P10 / P50 / P90

    dc_vec = [];
    if ~isfield(OUT,'deltap_cr_det') || isempty(OUT.deltap_cr_det), return; end
    dc_vec = OUT.deltap_cr_det(:);
    if ip == 0, return; end

    flds = {'mc_idx_fault_P10','mc_idx_fault_P50','mc_idx_fault_P90'};
    f = flds{ip};
    if ~isfield(OUT,f) || isempty(OUT.(f)), return; end
    i_star = OUT.(f)(wid, did);
    if ~isfinite(i_star), return; end
    if ~isfield(OUT,'deltap_cr_samples') || ~iscell(OUT.deltap_cr_samples), return; end

    nf  = numel(OUT.deltap_cr_samples);
    tmp = NaN(nf,1);
    for j = 1:nf
        s = OUT.deltap_cr_samples{j};
        if i_star <= numel(s), tmp(j) = s(i_star); end
    end
    dc_vec = tmp;
end


%% draw_heatmap_tab / draw_contour_tab
function draw_heatmap_tab(tab, data, d_list, well_list, well_labels, d_max, subtitle_str, COL, cbar_label, best_layout)
    if nargin < 9,  cbar_label  = '';  end
    if nargin < 10, best_layout = {}; end
    ax = axes('Parent', tab, 'Position', [0.10 0.12 0.75 0.75]); %#ok<LAXES>
    draw_heatmap(ax, data, d_list, well_list, well_labels, d_max, subtitle_str, COL, cbar_label, best_layout);
    add_unc_badge(tab, subtitle_str, COL);
end

function draw_contour_tab(tab, data, d_list, well_list, well_labels, d_max, subtitle_str, ...
                          cbar_label, round_levels, COL)
    ax = axes('Parent', tab, 'Position', [0.10 0.12 0.72 0.75]); %#ok<LAXES>
    draw_contour(ax, data, d_list, well_list, well_labels, d_max, subtitle_str, ...
        cbar_label, round_levels, COL);
    add_unc_badge(tab, subtitle_str, COL);
end

function add_unc_badge(parent, subtitle_str, COL)
    annotation(parent, 'textbox', [0.01 0.90 0.98 0.09], ...
        'String', subtitle_str, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 8.5, 'Color', [0.1 0.2 0.45], ...
        'BackgroundColor', [0.93 0.95 1.0], ...
        'EdgeColor', COL.accent, 'LineWidth', 0.8, ...
        'FitBoxToText', 'off', 'Interpreter', 'none');
end


%% draw_heatmap
function draw_heatmap(ax, data, d_list, well_list, well_labels, d_max, subtitle_str, COL, cbar_label, best_layout) %#ok<INUSL>
    if nargin < 9,  cbar_label  = '';  end
    if nargin < 10, best_layout = {}; end

    nw = length(well_list);
    nd = length(d_list);

    mask = false(nw, nd);
    for wi = 1:nw
        for di = 1:nd
            if isnan(data(wi,di))
                mask(wi,di) = true;
            end
        end
    end

    data_v = data;  data_v(mask) = NaN;
    dv_min = min(data_v(:), [], 'omitnan');
    dv_max = max(data_v(:), [], 'omitnan');
    if isnan(dv_min), return; end
    if dv_max - dv_min < 1e-10, dv_max = dv_min + 1; end

    cmap = make_blue_red_cmap(256);
    hold(ax, 'on');

    for wi = 1:nw
        for di = 1:nd
            val = data(wi, di);
            if mask(wi,di) || isnan(val)
                fill(ax, quad_x(di), quad_y(wi), ...
                    [0.87 0.87 0.89], 'EdgeColor','w','LineWidth',0.5);
                text(ax, di, wi, '—', 'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle','FontSize',8,'Color',[0.60 0.60 0.62]);
            else
                t   = min(1, max(0, (val-dv_min)/(dv_max-dv_min)));
                ci  = max(1, round(t*255)+1);
                clr = cmap(ci,:);
                fill(ax, quad_x(di), quad_y(wi), clr, 'EdgeColor','w','LineWidth',0.8);
                lum = 0.299*clr(1)+0.587*clr(2)+0.114*clr(3);
                text(ax, di, wi, sprintf('%.2f', val), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontSize',8,'FontWeight','bold','Color',ternary_color(lum));
            end
        end
    end

    for wi = 1:nw
        last_valid = find(~isnan(data(wi,:)), 1, 'last');
        if ~isempty(last_valid) && last_valid < nd
            plot(ax, [last_valid+0.5 last_valid+0.5], [wi-0.5 wi+0.5], '-', ...
                'Color', COL.warm, 'LineWidth', 3.2);
        end
    end

    colormap(ax, cmap);  caxis(ax, [dv_min dv_max]);
    xticks(ax, 1:nd);
    xticklabels(ax, arrayfun(@(v) sprintf('%.1f',v), d_list,'UniformOutput',false));
    yticks(ax, 1:nw);
    yticklabels(ax, well_labels);
    xlim(ax, [0.5 nd+0.5]);  ylim(ax, [0.5 nw+0.5]);
    xlabel(ax, 'Inter-well Distance  (km)', 'FontSize', 10);
    ylabel(ax, 'Number of Wells',            'FontSize', 10);
    title(ax,  'Heatmap table',             'FontSize', 11, 'FontWeight', 'bold');
    subtitle(ax, subtitle_str,              'FontSize', 8);
    cb = colorbar(ax);  cb.FontSize = 8;
    if ~isempty(cbar_label)
        cb.Label.String = cbar_label;
        cb.Label.FontSize = 9;
    end
    h_inf  = patch(ax,NaN,NaN,[0.87 0.87 0.89],'DisplayName','Infeasible','EdgeColor','none');
    h_dmax = plot(ax,NaN,NaN,'-','Color',COL.warm,'LineWidth',2.5,'DisplayName','d_{max}');
    legend(ax,[h_inf,h_dmax],'Location','southwest','FontSize',7);
    set(ax,'Box','on');
    hold(ax, 'off');

    setappdata(ax, 'hm_d_list',      d_list);
    setappdata(ax, 'hm_labels',      well_labels);
    setappdata(ax, 'hm_data',        data);
    setappdata(ax, 'hm_best_layout', best_layout);
    fig = ancestor(ax, 'figure');
    dcm = datacursormode(fig);
    dcm.UpdateFcn = @heatmap_tip_cb;
end


%% draw_contour
function draw_contour(ax, data, d_list, well_list, well_labels, d_max, subtitle_str, ...
                      cbar_label, round_levels, COL) %#ok<INUSL>

    nw = length(well_list);
    nd = length(d_list);
    data_plot = data;
    if ~any(isfinite(data_plot(:))), return; end

    finite_mask = isfinite(data_plot);
    if sum(finite_mask(:)) < 1
        text(ax, 0.5, 0.5, 'No finite data to display.', ...
            'Units','normalized','HorizontalAlignment','center', ...
            'FontSize', 9, 'Color', [0.45 0.45 0.45]);
        axis(ax, 'off');
        return;
    end

    hold(ax, 'on');

    h_img = imagesc(ax, [1, nd], [1, nw], data_plot);
    set(ax, 'YDir', 'normal');
    set(h_img, 'AlphaData', double(finite_mask));

    dv_min = min(data_plot(:), [], 'omitnan');
    dv_max = max(data_plot(:), [], 'omitnan');
    if dv_max - dv_min < 1e-10, dv_max = dv_min + 1; end

    colormap(ax, make_blue_red_cmap(256));
    caxis(ax, [dv_min dv_max]);
    cb = colorbar(ax);  cb.Label.String = cbar_label;  cb.FontSize = 8;

    y_idx = 1:nw;
    x_idx = 1:nd;

    finite_cols = any(finite_mask, 1);
    finite_rows = any(finite_mask, 2);
    if sum(finite_cols) >= 2 && sum(finite_rows) >= 2
        try
            [C, h] = contour(ax, x_idx, y_idx, data_plot, 30, ...
                'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);
            if round_levels, h.LevelList = round(h.LevelList, 2); end
            clabel(C, h, 'FontSize', 9, 'FontWeight', 'bold', ...
                'Color', [0.20 0.20 0.20], 'LabelSpacing', 400);
        catch
        end
    end

    d_boundary_idx = NaN(1, nw);
    for wi = 1:nw
        last_valid = find(~isnan(data(wi,:)), 1, 'last');
        if ~isempty(last_valid) && last_valid < nd
            d_boundary_idx(wi) = last_valid;
        end
    end
    valid_b = ~isnan(d_boundary_idx);
    if any(valid_b)
        plot(ax, d_boundary_idx(valid_b), y_idx(valid_b), '-', ...
            'Color', COL.warm, 'LineWidth', 2.4, 'DisplayName', 'd_{max}');
    end

    xlim(ax, [0.5, nd+0.5]);
    ylim(ax, [0.5, nw+0.5]);
    xticks(ax, x_idx);
    xticklabels(ax, arrayfun(@(v) sprintf('%.1f',v), d_list,'UniformOutput',false));
    yticks(ax, y_idx);
    yticklabels(ax, well_labels);
    xlabel(ax, 'Inter-well Distance  (km)', 'FontSize', 10);
    ylabel(ax, 'Number of Wells',            'FontSize', 10);
    title(ax,  'Contour plot',              'FontSize', 11, 'FontWeight', 'bold');
    subtitle(ax, subtitle_str,              'FontSize', 8);
    legend(ax, 'Location', 'best', 'FontSize', 7);
    set(ax, 'Color', COL.bg, 'Box', 'on', 'FontSize', 10);
    hold(ax, 'off');
end


%% fig_run_summary
function fig_run_summary(OUT, R, parent, COL)

    ax = axes('Parent', parent, 'Position', [0 0 1 1]); %#ok<LAXES>
    axis(ax, 'off');
    set(ax, 'Color', COL.bg);

    function s = yn(flag)
        if ischar(flag) || isstring(flag)
            if strcmpi(char(flag),'yes'), s = 'Yes';  else, s = 'No'; end
        elseif islogical(flag) || isnumeric(flag)
            if flag, s = 'Yes'; else, s = 'No'; end
        else
            s = '—';
        end
    end

    ref_pt   = char(getfield_safe(OUT, 'center_point_type', '—'));
    ft_raw   = char(getfield_safe(OUT, 'fault_type',        '—'));
    if strcmpi(ref_pt,'Fault')
        ft = ft_raw;
    else
        ft = 'N/A';
    end
    n_faults = 0;
    if isfield(OUT,'GEOM') && isfield(OUT.GEOM,'nr_fault')
        n_faults = OUT.GEOM.nr_fault;
    end

    stress_unc  = char(getfield_safe(OUT, 'stress_uncertain',  'no'));
    fault_unc   = char(getfield_safe(OUT, 'fault_uncertainty', 'no'));
    perm_unc    = OUT.perm_uncertain;
    n_MC_raw    = getfield_safe(OUT, 'n_MC', NaN);
    n_MC_set    = isfield(OUT,'n_MC_set') && OUT.n_MC_set;

    time_yr     = NaN;
    if isfield(OUT,'GEOM') && isfield(OUT.GEOM,'time')
        time_yr = OUT.GEOM.time / (86400*365);
    end
    correction  = '—';
    if isfield(OUT,'GEOM') && isfield(OUT.GEOM,'correction')
        correction = char(OUT.GEOM.correction);
    end

    use_sf_val = '—';
    if isfield(OUT,'GEOM') && isfield(OUT.GEOM,'use_sf')
        if OUT.GEOM.use_sf
            use_sf_val = 'Yes';
        else
            use_sf_val = 'No';
        end
    end

    r_inf_val = '—';
    if isfield(OUT,'GEOM') && isfield(OUT.GEOM,'R_influence')
        r_inf_km = OUT.GEOM.R_influence / 1000;
        r_inf_val = sprintf('%.1f  km', r_inf_km);
    end

    area_res  = getfield_safe(R, 'area_res', NaN);
    thick     = getfield_safe(R, 'thick',    NaN);
    por       = getfield_safe(R, 'por',      NaN);
    perm_ref  = getfield_safe(R, 'perm',     NaN);
    perm_mD   = perm_ref / 1e-15;

    btype = '—';
    if isfield(R, 'boundary') && ~isempty(R.boundary)
        btype = char(R.boundary);
    elseif isfield(OUT,'GEOM') && isfield(OUT.GEOM,'correction')
        btype = 'Closed (finite area)';
    end

    any_unc_active = strcmpi(stress_unc,'yes') || strcmpi(fault_unc,'yes') || ...
                     (islogical(perm_unc) && perm_unc) || ...
                     (isnumeric(perm_unc) && perm_unc > 0) || ...
                     (ischar(perm_unc) && strcmpi(perm_unc,'yes'));

    if any_unc_active
        if n_MC_set && ~isnan(n_MC_raw)
            n_MC_disp = fmt_val(n_MC_raw, '');
        else
            n_MC_disp = 'N/A';
        end
        agg_desc = 'Joint MC  —  P10 / P50 / P90 of Q_max across realizations';
    else
        n_MC_disp = 'N/A  —  no uncertainty sources active';
        agg_desc  = 'N/A  —  no uncertainty sources active';
    end

    rows = {
        '─── REFERENCE & FAULT ───────────────────────────────────────', '';
        'Reference point type',            ref_pt;
        'Fault type',                      ft;
        'Number of faults',                num2str(n_faults);
        '', '';
        '─── RESERVOIR GEOMETRY ──────────────────────────────────────', '';
        'Reservoir area',                  fmt_val(area_res, ' km²');
        'Reservoir thickness',             fmt_val(thick,    ' m');
        'Porosity',                        fmt_val(por,      '');
        'Reference permeability',          fmt_val(perm_mD,  ' mD');
        'Boundary condition',              btype;
        'Well correction factor',          correction;
        'Shape factor',                    use_sf_val;
        '', '';
        '─── UNCERTAINTY SOURCES ─────────────────────────────────────', '';
        'Stress uncertainty (stress MC)',  yn(stress_unc);
        'Fault geometry uncertainty',      yn(fault_unc);
        'Permeability uncertainty (MC)',   yn(perm_unc);
        'Monte Carlo realisations (n_MC)', n_MC_disp;
        'Uncertainty method',              agg_desc;
        '', '';
        '─── INJECTION SETTINGS ──────────────────────────────────────', '';
        'Injection period',                fmt_val(time_yr, ' years');
        'Radius of influence  R',          r_inf_val;
        'Well scenarios (n_wells)',         num2str(numel(OUT.well_list));
        'Distance scenarios (n_dist)',      num2str(numel(OUT.d_list));
        'Maximum number of wells',          num2str(max(OUT.well_list));
        'Max. sustainable inj. rate/well',  compute_max_qm_per_well(OUT);
    };

    text(ax, 0.50, 0.97, 'CO2BLOCK2  —  Run Summary', ...
        'Units','normalized','HorizontalAlignment','center', ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', COL.accent);
    text(ax, 0.50, 0.935, 'Audit trail of all active settings and uncertainty sources', ...
        'Units','normalized','HorizontalAlignment','center', ...
        'FontSize', 10, 'Color', COL.muted);

    x_lbl = 0.08;   x_val = 0.52;
    y0    = 0.90;    dy    = 0.030;

    for r = 1:size(rows, 1)
        lbl_r = rows{r,1};
        val_r = rows{r,2};
        y     = y0 - (r-1)*dy;

        if isempty(lbl_r) && isempty(val_r)
            continue;
        end

        if isempty(val_r)
            text(ax, x_lbl, y, lbl_r, 'Units','normalized', ...
                'FontSize', 8.5, 'FontWeight', 'bold', 'Color', COL.accent, ...
                'Interpreter','none');
        else
            text(ax, x_lbl, y, lbl_r, 'Units','normalized', ...
                'FontSize', 9, 'Color', [0.15 0.15 0.15], 'Interpreter','none');
            text(ax, x_val, y, val_r, 'Units','normalized', ...
                'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.05 0.25 0.50], ...
                'Interpreter','none');
        end
    end

    text(ax, 0.50, 0.02, ...
        'Settings sourced from General_settings.xlsx and Storage_unit_properties.xlsx', ...
        'Units','normalized','HorizontalAlignment','center', ...
        'FontSize', 8, 'Color', COL.muted, 'Interpreter','none');
end


function s = compute_max_qm_per_well(OUT)
    if isfield(OUT,'maxQ') && isnumeric(OUT.maxQ) && isfinite(OUT.maxQ) && OUT.maxQ > 0
        s = sprintf('%.4g  Mt/yr', OUT.maxQ);
        return;
    end

    try
        raw    = readcell('General_settings.xlsx', 'Sheet', 'General_settings', 'Range', 'A:B');
        raw    = raw(cellfun(@(x) ischar(x) && ~isempty(strtrim(x)), raw(:,1)), :);
        params = containers.Map(strtrim(raw(:,1)), raw(:,2));
        if isKey(params, 'maxQ')
            v = double(params('maxQ'));
            if isfinite(v) && v > 0
                s = sprintf('%.4g  Mt/yr', v);
                return;
            end
        end
    catch
    end

    s = '—';
end


function v = getfield_safe(s, fname, default)
    if isfield(s, fname) && ~isempty(s.(fname))
        v = s.(fname);
    else
        v = default;
    end
end

function s = fmt_val(v, unit)
    if isnumeric(v) && isfinite(v)
        if v == round(v)
            s = [num2str(v, '%g') unit];
        else
            s = [num2str(v, '%.4g') unit];
        end
    else
        s = '—';
    end
end


%% local helpers

function crit_arr = get_crit_ratio(OUT, wid, did)
    try
        crit_arr = OUT.crit_ratio_det{wid}{did}(:);
    catch
        crit_arr = [];
    end
end


function clr = ratio_traffic_color(r, COL)
    if isnan(r) || r <= 0
        clr = [0.80 0.80 0.82];
    elseif r >= 1.0
        clr = COL.warm;
    elseif r >= 0.75
        clr = COL.amber;
    else
        clr = COL.mid;
    end
end

function tc = ternary_color(lum)
    if lum < 0.52, tc = [1 1 1]; else, tc = [0.08 0.08 0.08]; end
end

function xs = quad_x(di)
    xs = [di-0.5, di+0.5, di+0.5, di-0.5];
end

function ys = quad_y(row)
    ys = [row-0.5, row-0.5, row+0.5, row+0.5];
end

function cmap = make_fault_cmap(n)
    lo=[0.13 0.70 0.20]; mid=[0.97 0.73 0.00]; hi=[0.84 0.10 0.10];
    h1=floor(n/2); h2=n-h1;
    cmap=[interp_ramp(lo,mid,h1); interp_ramp(mid,hi,h2)];
end

function cmap = make_blue_red_cmap(n)
    lo=[0.17 0.44 0.70]; mid=[0.97 0.97 0.97]; hi=[0.84 0.19 0.15];
    h1=floor(n/2); h2=n-h1;
    cmap=[interp_ramp(lo,mid,h1); interp_ramp(mid,hi,h2)];
end

function cmap = make_pressure_cmap(n)
% Dark blue (high ΔP) → steel blue → near-white (low ΔP)
    lo  = [0.03 0.07 0.35];   % dark navy
    mid = [0.13 0.53 0.80];   % steel blue
    hi  = [0.92 0.96 1.00];   % near-white
    h1 = floor(n/2); h2 = n - h1;
    cmap = [interp_ramp(hi, mid, h1); interp_ramp(mid, lo, h2)];
end

function c = interp_ramp(a, b, n)
    c=[linspace(a(1),b(1),n)', linspace(a(2),b(2),n)', linspace(a(3),b(3),n)'];
end


%% compute_pressure_map / pressure_at_points
function [p_map, X, Y] = compute_pressure_map(GEOM, R, wid, did, Q_Mtyr)
% Overpressure map [MPa].  Q_Mtyr = per-well rate [Mt/yr]; omit or leave
% empty for the natural reference rate Q0 = M0/w (original behaviour).

    if nargin < 5, Q_Mtyr = []; end
    grid_div = 200;
    [X, Y] = meshgrid(linspace(0, R.x_length, grid_div), ...
                      linspace(0, R.y_length, grid_div));
    p_map = reshape(pressure_at_points(GEOM, R, wid, did, Q_Mtyr, X(:), Y(:)), size(X));
end


function p_pts = pressure_at_points(GEOM, R, wid, did, Q_Mtyr, px, py)
% Superposed Nordbotten overpressure [MPa] at arbitrary points.
%
% Pressure is NOT proportional to Q: the rate sets both the amplitude p_c
% and the plume radius through csi / psi.  Every quantity that depends on Q
% is therefore recomputed here from the requested rate.

    perm=R.perm; visc_w=R.visc_w; compr=R.compr; rc=R.rc;
    gamma=R.gamma; delta=R.delta; omega=R.omega;
    dens_c=R.dens_c; thick=R.thick; por=R.por;
    area_res=R.area_res;

    rw=GEOM.rw; time=GEOM.time; correction=GEOM.correction;
    use_sf = isfield(GEOM,'use_sf') && GEOM.use_sf;

    w        = GEOM.well_list(wid);
    distance = GEOM.distance_matrix(wid,did);
    wx = GEOM.well_coords_x{wid,did}(:);
    wy = GEOM.well_coords_y{wid,did}(:);

    R_influence = sqrt(2.246*perm*time/(visc_w*compr));
    M0 = perm/1e-13;

    % per-well rate [m3/s]
    if isempty(Q_Mtyr)
        Q = M0*1e9/dens_c/365/86400/w;      % natural reference rate
    else
        Q = Q_Mtyr*1e9/dens_c/365/86400;    % requested rate
    end

    csi = sqrt(Q*time/pi/por/thick);
    psi = exp(omega)*csi;
    p_c = (Q*visc_w)/(2*pi*thick*perm)/1e6;

    sup_err = 0;
    if strcmp(correction,'on') && w >= 9 && R_influence*csi/distance^2 >= 1
        sup_err = w*delta/4*log(R_influence*csi/distance^2);
    end

    if use_sf
        C_A_v = GEOM.shape_CA_cache{wid,did};
        Rx_v  = GEOM.shape_Rext_cache{wid,did};
    else
        C_A_v = [];  Rx_v = [];
    end

    px = px(:);  py = py(:);
    p_pts = zeros(numel(px),1);
    for k = 1:numel(px)
        r = sqrt((wx-px(k)).^2 + (wy-py(k)).^2);
        r(r==0) = rw;
        PD = Nordbotten_solution(r, R_influence, psi, rc, gamma, ...
                                 use_sf, C_A_v, Rx_v, area_res);
        p_pts(k) = sum(PD)*p_c - sup_err*p_c;
    end
end


%% heatmap_tip_cb
function txt = heatmap_tip_cb(~, evt)
    ax = ancestor(evt.Target, 'axes');
    if isempty(ax) || ~isappdata(ax, 'hm_d_list')
        txt = '';
        return;
    end
    d_list      = getappdata(ax, 'hm_d_list');
    labels      = getappdata(ax, 'hm_labels');
    data        = getappdata(ax, 'hm_data');
    best_layout = getappdata(ax, 'hm_best_layout');
    pos = evt.Position;
    di = round(pos(1));
    wi = round(pos(2));
    if di < 1 || di > numel(d_list) || wi < 1 || wi > numel(labels)
        txt = '';
        return;
    end
    val = data(wi, di);

    grid_str = '';
    if ~isempty(best_layout) && wi <= size(best_layout,1) && di <= size(best_layout,2)
        grid_str = best_layout{wi, di};
    end

    if isnan(val)
        txt = {sprintf('Wells: %s', labels{wi}), ...
               sprintf('d = %.1f km', d_list(di)), ...
               'Infeasible'};
    else
        if ~isempty(grid_str)
            txt = {sprintf('Wells: %s', labels{wi}), ...
                   sprintf('d = %.1f km', d_list(di)), ...
                   sprintf('Grid: %s', grid_str), ...
                   sprintf('Value: %.4f', val)};
        else
            txt = {sprintf('Wells: %s', labels{wi}), ...
                   sprintf('d = %.1f km', d_list(di)), ...
                   sprintf('Value: %.4f', val)};
        end
    end
end