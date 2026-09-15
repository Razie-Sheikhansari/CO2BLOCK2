function [OUT, R] = calculate(S, R, G)
% Runs the full CO2BLOCK pipeline in sequence.
% Fluid props -> MC sampling -> critical overpressure -> pressure buildup
% -> fault reactivation -> storage capacity (det + fault det + P10/P50/P90).

 
    R = compute_fluid_properties(R, G);
    MC = run_monte_carlo(S, R, G);
    GEO = compute_critical_overpressure(S, G, R, MC);
    [P, GEOM] = compute_pressure_buildup(S, R, GEO);
    P.perm_ref = R.perm;
    % fault reactivation criticality ratio
    CR = compute_reactivation_risk(S, P, GEO, MC);
    % storage capacity - det, fault-constrained det, joint MC percentiles
    OUT = compute_storage_capacity(S, R, P, CR, MC, GEO, GEOM);
end


% 5 scenarios
% Det      : Q_max at ref perm, mean stress, no fault constraint
% Fault Det: Q_max at ref perm, mean stress, most critical fault
% Joint P10: 10th percentile of Q_max across MC realizations
% Joint P50: 50th percentile of Q_max across MC realizations
% Joint P90: 90th percentile of Q_max across MC realizations

