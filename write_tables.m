function write_tables(OUT, S)
% Write CO2BLOCK2 injection rate and storage capacity tables to .xls files.
% Desired mode writes two file pairs: one at the user's (x,y) point, one at the central reference.

    is_fault   = OUT.center_point_type == "Fault";
    is_desired = OUT.center_point_type == "Desired";
    any_mc     = isfield(OUT,'Q_M_each_joint_P50') && ~isempty(OUT.Q_M_each_joint_P50);

    d_list      = OUT.d_list;
    x_grid_list = OUT.x_grid_list;
    y_grid_list = OUT.y_grid_list;

    fprintf('\n=== Writing output tables (write_tables.m v5 - Desired mode fix) ===\n');
    fprintf('  center_point_type : %s\n', OUT.center_point_type);
    fprintf('  Available OUT fields: %s\n', strjoin(fieldnames(OUT), ', '));

    %% desired-point result (Desired mode only)
    % Q_M_each_det / V_M_det hold the result at the user's (x,y) coordinates
    if is_desired
        write_pair(OUT, {'Q_M_each_det'}, {'V_M_det'}, ...
            d_list, x_grid_list, y_grid_list, ...
            'Q_M_desired_max_per_well_inj_rate.xls', ...
            'V_M_desired_max_storage_capacity.xls');
    end

    %% deterministic central reference (always)
    % Q_M_each_ref / V_M_ref use the geometric-centre pressure field regardless of mode
    % falls back to old field names (_2) for backward compatibility
    write_pair(OUT, {'Q_M_each_ref','Q_M_each_2'}, ...
                    {'V_M_ref',     'V_M_2'}, ...
        d_list, x_grid_list, y_grid_list, ...
        'Q_M_det_max_per_well_inj_rate.xls', ...
        'V_M_det_max_storage_capacity.xls');

    %% fault-constrained deterministic (Fault mode only)
    if is_fault
        write_pair(OUT, {'Q_M_each_fault_det','Q_M_each_fault_2'}, ...
                        {'V_M_fault_det',     'V_M_fault_2'}, ...
            d_list, x_grid_list, y_grid_list, ...
            'Q_M_fault_det_max_per_well_inj_rate.xls', ...
            'V_M_fault_det_max_storage_capacity.xls');
    end

    %% joint MC P10 / P50 / P90
    if any_mc
        % reservoir-limited joint
        if ~is_fault
            for pct = {'P10','P50','P90'}
                p = pct{1};
                write_pair(OUT, {['Q_M_each_joint_' p]}, ...
                                {['V_M_joint_' p]}, ...
                    d_list, x_grid_list, y_grid_list, ...
                    ['Q_M_joint_' p '_max_per_well_inj_rate.xls'], ...
                    ['V_M_joint_' p '_max_storage_capacity.xls']);
            end
        end
        % fault-constrained joint (Fault mode only)
        if is_fault
            for pct = {'P10','P50','P90'}
                p = pct{1};
                write_pair(OUT, {['Q_M_each_joint_fault_' p]}, ...
                                {['V_M_joint_fault_' p]}, ...
                    d_list, x_grid_list, y_grid_list, ...
                    ['Q_M_joint_fault_' p '_max_per_well_inj_rate.xls'], ...
                    ['V_M_joint_fault_' p '_max_storage_capacity.xls']);
            end
        end
    end

    fprintf('=== Done ===\n');
end


function write_pair(OUT, qfields, vfields, d_list, x_grid_list, y_grid_list, qfile, vfile)
% Build Q and V tables from raw [nw x nd] matrices and write to .xls files.
% qfields / vfields: cell arrays of candidate field names, tried in order.
    Q = safe_get(OUT, qfields);
    V = safe_get(OUT, vfields);
    % Lambert-W / Nordbotten can produce complex Q_M values near the branch cut;
    % take the real part (same as CO2BLOCK_plots)
    if ~isempty(Q), Q = real(Q); end
    if ~isempty(V), V = real(V); end
    if isempty(Q) || isempty(V)
        if iscell(qfields), qn = strjoin(qfields,', '); else, qn = qfields; end
        if iscell(vfields), vn = strjoin(vfields,', '); else, vn = vfields; end
        fprintf('  Skipped (no data): {%s} / {%s}\n', qn, vn);
        return;
    end

    % diagnostics: report finite / NaN / complex counts
    n_total     = numel(Q);
    n_complex_Q = nnz(imag(Q(:)) ~= 0);
    n_complex_V = nnz(imag(V(:)) ~= 0);
    Q = real(Q);  % belt-and-suspenders (also done above)
    V = real(V);
    n_finite_Q = nnz(isfinite(Q(:)));
    n_finite_V = nnz(isfinite(V(:)));

    fprintf('  [%s]  %d/%d finite values', qfile, n_finite_Q, n_total);
    if n_complex_Q > 0, fprintf(' (%d were complex - real part kept)', n_complex_Q); end
    fprintf('\n');
    fprintf('  [%s]  %d/%d finite values', vfile, n_finite_V, n_total);
    if n_complex_V > 0, fprintf(' (%d were complex - real part kept)', n_complex_V); end
    fprintf('\n');

    if n_finite_Q == 0
        fprintf('  *** WARNING: %s has NO finite data - all cells are NaN.\n', qfile);
        fprintf('      This means the Lambert-W equation has no real solution\n');
        fprintf('      for ANY grid config at the deterministic overpressure.\n');
        fprintf('      The Excel file will have blank data cells.  Check if\n');
        fprintf('      the Joint-P10/P50/P90 Excel files have data instead.\n');
    end

    TQ = to_table(Q, d_list, x_grid_list, y_grid_list);
    TV = to_table(V, d_list, x_grid_list, y_grid_list);
    writetable(TQ, qfile, 'WriteRowNames', true);
    writetable(TV, vfile, 'WriteRowNames', true);
    fprintf('  Saved: %s\n', qfile);
    fprintf('  Saved: %s\n', vfile);
end

function T = to_table(M, d_list, x_grid_list, y_grid_list)
    vn = cellstr(arrayfun(@(d) sprintf('d_%.0f_m', d*1000), d_list, 'UniformOutput', false));
    rn = arrayfun(@(nx,ny) sprintf('%dx%d', nx, ny), ...
        x_grid_list(:), y_grid_list(:), 'UniformOutput', false);
    T  = array2table(M, 'VariableNames', vn, 'RowNames', rn);
    T.Properties.DimensionNames{1} = 'grid_layout';
end

function v = safe_get(OUT, fnames)
% Try each field name in turn; return the first non-empty match.
    if ischar(fnames), fnames = {fnames}; end
    for i = 1:numel(fnames)
        if isfield(OUT, fnames{i}) && ~isempty(OUT.(fnames{i}))
            fprintf('    safe_get: matched field ''%s'' [%dx%d]\n', ...
                fnames{i}, size(OUT.(fnames{i}),1), size(OUT.(fnames{i}),2));
            v = OUT.(fnames{i});
            return;
        else
            fprintf('    safe_get: tried ''%s'' - %s\n', fnames{i}, ...
                ternary(isfield(OUT,fnames{i}), 'exists but empty', 'not found'));
        end
    end
    fprintf('    safe_get: NO match found among {%s}\n', strjoin(fnames,', '));
    v = [];
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end