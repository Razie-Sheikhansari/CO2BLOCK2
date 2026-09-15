function [Sigma_n, Tau] = stress_projection(delta, SH_max, Sh_min, Sv, Azi, Dip)
% Project principal stresses onto a fault plane and return normal and total shear stress.
% delta: azimuth of SHmax from North [deg]. Azi, Dip: fault strike and dip [deg].

    % rotation matrix from principal stress frame (SHmax, Shmin, Sv) to geographic N-E-Down
    % SHmax lies at azimuth delta from North; Sv is vertical (assumed principal)
    L_PG = [cos(deg2rad(delta))  sin(deg2rad(delta))  0; ...
            sin(deg2rad(delta)) -cos(deg2rad(delta))  0; ...
            0                    0                   -1];

    % principal stress tensor (diagonal in principal frame)
    Sigma = [SH_max 0 0; 0 Sh_min 0; 0 0 Sv];

    % stress tensor in geographic N-E-Down frame
    Sigma_G = L_PG * Sigma * L_PG';

    % fault unit vectors in N-E-Down
    n_n = [-(sin(deg2rad(Azi))) * (sin(deg2rad(Dip))); ...   % fault normal (upward, toward footwall)
            (cos(deg2rad(Azi))) * (sin(deg2rad(Dip))); ...
           -cos(deg2rad(Dip))];

    n_s = [cos(deg2rad(Azi)); sin(deg2rad(Azi)); 0];          % along-strike direction

    n_d = [-(sin(deg2rad(Azi))) * (cos(deg2rad(Dip))); ...   % down-dip direction (third component + = downward)
            (cos(deg2rad(Azi))) * (cos(deg2rad(Dip))); ...
            sin(deg2rad(Dip))];

    % traction on fault plane, then resolve into normal and shear components
    Traction   = Sigma_G * n_n;                              % traction vector in N-E-Down
    Sigma_n    = transpose(n_n) * Traction;                  % normal stress on fault
    Tau_s      = transpose(n_s) * Traction;                  % shear along strike
    Tau_d      = transpose(n_d) * Traction;                  % shear along dip (+ = downward)
    Tau_vector = Tau_s * n_s + Tau_d * n_d;                 % total shear stress vector
    Tau        = norm(Tau_vector);                           % total shear stress magnitude
end