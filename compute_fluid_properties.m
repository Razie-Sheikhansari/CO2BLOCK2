function R = compute_fluid_properties(R, G)
% Fill reservoir fluid properties using eos.
% T_mean and pres_mean resolved from: known value > gradient > default.
% eos is only called here; direct inputs in xlsx override it.
%
%  Default values:
%     pres_grad   10.2    MPa/km   hydrostatic
%     temp_surf   15    °C
%     temp_grad   33    °C/km
%     depth       =     depth_mean − thick/2
%     depth_water  0    m  (onshore)
%     cr          5e-4  MPa⁻¹
%     cw          3e-4  MPa⁻¹
%     dens_c, visc_c, visc_w from eos(T_mean, pres_mean, salinity)

    if G.depth_mean == 0 || isnan(G.depth_mean)
        error(['compute_fluid_properties: depth_mean is not set in ' ...
               'Geomech_parameters.xlsx.\n' ...
               'This is a required input — no default exists for reservoir depth.']);
    end

    depth_water = G.depth_water;   % 0 = onshore

    % top reservoir depth — default to mid-depth minus half thickness
    if G.depth == 0 || isnan(G.depth)
        G.depth = G.depth_mean - R.thick / 2;
    end

    % temperature
    have_T_direct = isfield(G,'T_mean_known') && isfinite(G.T_mean_known);
    if have_T_direct
        T_mean = G.T_mean_known;
        fprintf('compute_fluid_properties: T_mean = %.1f °C  (from T_mean_known)\n', T_mean);
    else
        depth_eff_T = G.depth_mean - depth_water;
        T_mean      = G.temp_surf + (G.temp_grad / 1e3) * depth_eff_T;
        fprintf('compute_fluid_properties: T_mean = %.1f °C  (surf=%.1f + %.1f°C/km × %.0fm)\n', ...
            T_mean, G.temp_surf, G.temp_grad, depth_eff_T);
    end

    % pressure
    have_P_direct = isfield(G,'pres_mean_known') && isfinite(G.pres_mean_known);
    if have_P_direct
        pres_mean = G.pres_mean_known;
        fprintf('compute_fluid_properties: pres_mean = %.2f MPa  (from pres_mean_known)\n', pres_mean);
    else
        depth_eff_P = G.depth_mean - depth_water;
        pres_mean   = G.pres_offset + G.pres_grad * depth_eff_P;
        fprintf('compute_fluid_properties: pres_mean = %.2f MPa  (offset=%.2f + %.1fMPa/km × %.0fm)\n', ...
            pres_mean, G.pres_offset, G.pres_grad * 1e3, depth_eff_P);
    end

    % fluid properties from eos (skipped if provided directly)
    if R.dens_c == 0 || isnan(R.dens_c)
        [~, R.dens_c, ~] = eos(T_mean, pres_mean, R.salinity, 0);
        fprintf('compute_fluid_properties: dens_c = %.1f kg/m³ (eos)\n', R.dens_c);
    else
        fprintf('compute_fluid_properties: dens_c = %.1f kg/m³ (input)\n', R.dens_c);
    end

    if R.visc_c == 0 || isnan(R.visc_c)
        [~, ~, R.visc_c] = eos(T_mean, pres_mean, R.salinity, R.dens_c);
        fprintf('compute_fluid_properties: visc_c = %.4e Pa.s (eos)\n', R.visc_c);
    else
        fprintf('compute_fluid_properties: visc_c = %.4e Pa.s (input)\n', R.visc_c);
    end

    if R.visc_w == 0 || isnan(R.visc_w)
        [R.visc_w, ~, ~] = eos(T_mean, pres_mean, R.salinity, 0);
        fprintf('compute_fluid_properties: visc_w = %.4e Pa.s (eos)\n', R.visc_w);
    else
        fprintf('compute_fluid_properties: visc_w = %.4e Pa.s (input)\n', R.visc_w);
    end

    % compressibility defaults (cr=5e-4, cw=3e-4 MPa⁻¹)
    if R.cr == 0 || isnan(R.cr)
        R.cr = 5e-4 / 1e6;
        fprintf('compute_fluid_properties: cr defaulted to 5e-4 MPa⁻¹\n');
    end
    if R.cw == 0 || isnan(R.cw)
        R.cw = 3e-4 / 1e6;
        fprintf('compute_fluid_properties: cw defaulted to 3e-4 MPa⁻¹\n');
    end

    % viscosity ratios and bulk compressibility for Nordbotten solution
    R.gamma = R.visc_c / R.visc_w;
    R.delta = (R.visc_w - R.visc_c) / R.visc_w;
    R.omega = (R.visc_c + R.visc_w) / (R.visc_c - R.visc_w) * ...
               log(sqrt(R.visc_c / R.visc_w)) - 1;
    R.compr = (R.cr + R.cw) * R.por;

    R.pres_mean  = pres_mean;
    R.T_mean     = T_mean;
    R.depth      = G.depth;
    R.depth_mean = G.depth_mean;

    fprintf('compute_fluid_properties: done | T=%.1f°C | P=%.2f MPa | dens_c=%.1f | visc_c=%.4e | visc_w=%.4e\n', ...
        T_mean, pres_mean, R.dens_c, R.visc_c, R.visc_w);
end