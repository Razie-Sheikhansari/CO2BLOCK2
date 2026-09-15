function PD = Nordbotten_solution(r, R, psi, R_ext, gamma, use_sf, C_A, R_ext_well, area_res)
% Dimensionless pressure at r using the Nordbotten analytical solution.
% Routes to the shape-factor variant when use_sf is true.

    if nargin < 6 || ~use_sf
        PD = nordbotten_standard(r, R, psi, R_ext, gamma);
    else
        PD = nordbotten_shape_factor(r, R, psi, R_ext, gamma, C_A, R_ext_well, area_res);
    end
end

%% standard path
function PD = nordbotten_standard(r, R, psi, R_ext, gamma)
    r  = max(r, 1.0);
    PD = zeros(size(r));

    if R >= psi
        % pressure front contains plume
        FD_psi    = FD_Nor(psi, R, R_ext);
        mask1     = r <= psi;
        PD(mask1) = gamma * log(psi ./ r(mask1)) + FD_psi;
        mask2     = ~mask1 & (r <= R);
        PD(mask2) = FD_Nor(r(mask2), R, R_ext);

    else
        % plume has grown past the pressure front (R < psi)
        % gamma = visc_c/visc_w < 1, so r_bot = psi*sqrt(gamma) < psi < r_top = psi/sqrt(gamma)
        % three zones: CO2 log rise (r <= r_bot), linear ramp, brine (r > r_top)
        r_bot     = psi * sqrt(gamma);
        r_top     = psi / sqrt(gamma);
        FD_rtop   = FD_Nor(r_top, R, R_ext);
        ramp_full = sqrt(gamma) / psi * (r_top - r_bot);

        mask1     = r <= r_bot;
        PD(mask1) = gamma * log(r_bot ./ r(mask1)) + ramp_full + FD_rtop;

        mask2     = r > r_bot & r <= r_top;
        PD(mask2) = sqrt(gamma) / psi * (r_top - r(mask2)) + FD_rtop;

        mask3     = r > r_top & r < R_ext;
        PD(mask3) = FD_Nor(r(mask3), R, R_ext);
    end
end

%% shape-factor path
function PD = nordbotten_shape_factor(r, R, psi, rc, gamma, C_A, R_ext_well, area_res)
    r  = max(r, 1.0);
    nw = numel(r);
    PD = zeros(nw, 1);

    bad_CA = ~isfinite(C_A) | C_A <= 0;
    if any(bad_CA)
        warning('Nordbotten_solution: %d/%d wells have invalid C_A - check x_length/y_length in Excel.', ...
            sum(bad_CA), nw);
    end

    % print key scalars once on first call
    persistent sf_called;
    if isempty(sf_called)
        fprintf('  Nordbotten_sf first call: R=%.1fm  psi=%.1fm  rc=%.1fm  R>rc=%d\n', ...
            R, psi, rc, R > rc);
        sf_called = true;
    end

    % r_bot and r_top are the same for all wells (R, psi, gamma are scalars)
    r_bot = psi * sqrt(gamma);
    r_top = psi / sqrt(gamma);

    for iw = 1:nw
        R_ext_iw = R_ext_well(iw);
        C_A_iw   = C_A(iw);
        r_iw     = r(iw);

        % fall back to standard FD if C_A is invalid
        if bad_CA(iw)
            fd = @(x) FD_Nor(x, R, rc);
        else
            fd = @(x) FD_Nor_sf(x, R, R_ext_iw, rc, area_res, C_A_iw);
        end

        if R >= psi
            % Case A: pressure front contains plume
            FD_psi_iw = fd(psi);
            if r_iw <= psi
                PD(iw) = gamma * log(psi / r_iw) + FD_psi_iw;
            elseif r_iw <= R && r_iw < R_ext_iw
                PD(iw) = fd(r_iw);
            end

        else
            % Case B: plume past pressure front
            FD_rtop   = fd(r_top);
            ramp_full = sqrt(gamma) / psi * (r_top - r_bot);

            if r_iw <= r_bot
                PD(iw) = gamma * log(r_bot / r_iw) + ramp_full + FD_rtop;
            elseif r_iw <= r_top
                PD(iw) = sqrt(gamma) / psi * (r_top - r_iw) + FD_rtop;
            elseif r_iw < R_ext_iw
                PD(iw) = fd(r_iw);
            end
        end
    end
end

%% FD_Nor_sf: shape-factor dimensionless pressure (scalar input)
% Returns 0 outside the pressure front. Branches on R vs rc, not R vs R_ext_well.
function FD = FD_Nor_sf(x, R, R_ext, R_ext_c, area_res, C_A)
    x  = max(x, 1.0);
    FD = 0;
    if x >= R || x >= R_ext   % outside pressure front or drainage area
        return;
    end
    if R <= R_ext_c          % radial log (same as standard)
        FD = log(R / x);
    else                     % shape-factor formula
        FD = (2*pi*R^2) / (2.25 * area_res * 1e6) + ...
             0.5 * log(4 * area_res * 1e6 / (x^2 * C_A * 1.781));
    end
end