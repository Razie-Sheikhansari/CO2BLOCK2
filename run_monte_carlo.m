function MC = run_monte_carlo(S, R, G)
% Draw MC samples for permeability and stress.
% Produces raw samples only; Layer 5 uses them directly without further aggregation.

    n_MC             = S.n_MC;
    stress_uncertain = S.stress_uncertain;

    MC.n_MC           = n_MC;
    MC.stress_active  = false;
    MC.stress_samples = struct();
    MC.perm_active    = false;
    MC.perm_samples   = [];

    % stress
    if stress_uncertain == "yes"
        depth_eff = G.depth_seismic - G.depth_water;
        pd_p     = build_dist(G.f_dist_p,     G.pres_offset,  depth_eff, ...
                              G.pres_grad_seismic,  G.a_p_grad/1000);
        pd_Sv    = build_dist(G.f_dist_Sv,    G.Sv_offset,    depth_eff, ...
                              G.Sv_grad_seismic,    G.a_Sv_grad/1000);
        pd_Shmin = build_dist(G.f_dist_Shmin, G.Shmin_offset, depth_eff, ...
                              G.Shmin_grad_seismic, G.a_Shmin_grad/1000);
        pd_SHmax = build_dist(G.f_dist_SHmax, G.SHmax_offset, depth_eff, ...
                              G.SHmax_grad_seismic, G.a_SHmax_grad/1000);
        pd_dir   = build_dist_angle(G.f_dist_dir, G.SHmax_dir_seismic, G.a_SHmax_dir);

        ss.r_p     = random(pd_p,     n_MC, 1);
        ss.r_Sv    = random(pd_Sv,    n_MC, 1);
        ss.r_Shmin = random(pd_Shmin, n_MC, 1);
        ss.r_SHmax = random(pd_SHmax, n_MC, 1);
        ss.r_dir   = random(pd_dir,   n_MC, 1);

        MC.stress_active  = true;
        MC.stress_samples = ss;
        fprintf('  stress MC  (%d draws)\n', n_MC);
    end

    % permeability
    if S.perm_uncertain == "yes"
        perm_mD  = R.perm / 1e-15;        % m2 -> mD
        log_perm = log10(perm_mD);         % sample in log10(mD) space
        switch upper(R.perm_dist)
            case 'NOR'
                pd = makedist('Normal', 'mu', log_perm, 'sigma', R.perm_std);
            case 'UNI'
                pd = makedist('Uniform', ...
                    'lower', log_perm - R.perm_std, ...
                    'upper', log_perm + R.perm_std);
            otherwise
                error('run_monte_carlo: unknown perm_dist "%s". Use "Nor" or "Uni".', ...
                    R.perm_dist);
        end
        perm_samples = (10.^random(pd, n_MC, 1)) * 1e-15;  % log10(mD) -> m2
        MC.perm_active  = true;
        MC.perm_samples = perm_samples;
        fprintf('  perm MC  (%d draws)  ref=%.4g mD  P10=%.4g  P90=%.4g\n', ...
            n_MC, R.perm/1e-15, prctile(perm_samples,10)/1e-15, prctile(perm_samples,90)/1e-15);
    end
end


function pd = build_dist(dist_type, offset, depth, grad, a_grad)
% Uniform or normal distribution centred at offset + depth*grad.
% half-width = depth * a_grad (uncertainty scales with depth).
    mean_val = offset + depth * grad;
    half     = depth * a_grad;
    if strcmp(dist_type, 'Uni')
        pd = makedist('Uniform', 'lower', mean_val-half, 'upper', mean_val+half);
    else
        pd = makedist('Normal', 'mu', mean_val, 'sigma', half);
    end
end


function pd = build_dist_angle(dist_type, center, half_range)
% Uniform or normal distribution for SHmax azimuth.
    if strcmp(dist_type, 'Uni')
        pd = makedist('Uniform', 'lower', center-half_range, 'upper', center+half_range);
    else
        pd = makedist('Normal', 'mu', center, 'sigma', half_range);
    end
end