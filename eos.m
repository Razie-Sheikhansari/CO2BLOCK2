% Brine viscosity from Batzle & Wang (1992).
% CO2 density from Redlich-Kwong (1949) with parameters from Spycher et al. (2003).
% CO2 viscosity from Altunin & Sakhabetdinov (1972).

function [brineviscosity, co2density, co2viscosity] = eos(T, p, salinity, co2dens)
% T [°C], p [MPa], salinity [ppm/1e6]

% brine viscosity
    brineviscosity = (0.1 + 0.333*salinity + (1.65 + 91.9*salinity^3) ...
        * exp(-(0.42*(salinity^0.8 - 0.17)^2 + 0.045)*T^0.8)) / 1e3;  % [Pa s]

    T = T + 273.15;   % °C → K
    p = p * 1e6;      % MPa → Pa

% CO2 density (RK EOS)
    a0 = 7.54;       % [Pa m6 K^0.5 mol^-2]
    a1 = -4.13e-3;   % [Pa m6 K^-0.5 mol^-2]
    b  = 2.78e-5;    % [m3 mol^-1]
    a  = a0 + a1*T;
    R  = 8.314472;   % [m3 Pa K^-1 mol^-1]

    A = -(R*T/p);
    B = -(R*T*b/p - a/(p*sqrt(T)) + b^2);
    C = -(a*b/(p*sqrt(T)));

    V_all  = roots([1 A B C]);
    V_real = V_all(imag(V_all) == 0 & real(V_all) > 0);
    V = max(V_real);            % supercritical CO2: one real root
    co2density = 0.044 / V;    % [kg/m3]

% CO2 viscosity (Altunin & Sakhabetdinov 1972)
    a10 =  0.248566120;
    a11 =  0.004894942;
    a20 = -0.373300660;
    a21 =  1.22753488;
    a30 =  0.363854523;
    a31 = -0.774229021;
    a40 = -0.0639070755;
    a41 =  0.142507049;

    Tr     = T / 304;        % reduced T  (Tc = 304 K)
    dens_r = co2dens / 468;  % reduced ρ  (ρc = 468 kg/m3)

    mu_0 = Tr^0.5 * (27.2246461 - 16.6346068/Tr + 4.66920556/Tr^2) * 1e-6;  % [Pa s]

    co2viscosity = double(mu_0 * exp( ...
        a10*dens_r   + a11*dens_r/Tr   + ...
        a20*dens_r^2 + a21*dens_r^2/Tr + ...
        a30*dens_r^3 + a31*dens_r^3/Tr + ...
        a40*dens_r^4 + a41*dens_r^4/Tr));  % [Pa s]

end