% CO2BLOCK2.m  —  Main script
%
% Estimates CO2 storage capacity of a saline reservoir for different
% combinations of injection well number and inter-well spacing, accounting
% for pressure buildup and fault reactivation risk

clearvars; close all;

%% 1. READ INPUTS
S = read_settings();       % General_settings.xlsx      
R = read_storage_unit();   % Storage_unit_properties.xlsx >> reservoir properties (raw)
G = read_geomech();        % Geomech_parameters.xlsx + Faults.xlsx

%% 2. RUN CALCULATION
%   calculate returns the completed R (with compr, gamma, delta, omega etc.
%   filled in by compute_fluid_properties) 
[OUT, R] = calculate(S, R, G);

%% 3. SAVE RESULTS
save('CO2BLOCK2_results.mat', 'OUT', 'S', 'R', 'G');

%% 4. WRITE OUTPUT TABLES
write_tables(OUT, S);

%% 5. PLOTS
CO2BLOCK2_plots(OUT, R);

fprintf('\nCO2BLOCK2 analysis complete.\n');
