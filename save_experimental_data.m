% save_data.m
Xexp = [0.334, 0.323, 0.337, 0.350, 0.365]; % Row vector
T = 273.15 + [47.7, 43.7, 48.4, 52.4, 57.4]; % Convert to Kelvin
Hfus = 25.4 * 1000; % [J/mol]
Tfus = 158.7 + 273.15; % [K]
R = 8.314; % Gas constant in J/(mol*K)
v1 = 85.54;
v2 = 96;

% Save the data to a .mat file
save('data.mat', 'Xexp', 'T', 'Hfus', 'Tfus', 'R', 'v1', 'v2');
