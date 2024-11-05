% Experimental Data Points for Solute 1
Xexp = [0.334, 0.323, 0.337, 0.350, 0.365]; % Row vector
T = 273.15 + [47.7, 43.7, 48.4, 52.4, 57.4]; % Convert to Kelvin

Hfus = 25.4 * 1000; % [J/mol]
Tm = 158.7 + 273.15; % [K] Melting point of the solute
R = 8.314; % Gas constant in J/(mol*K)

% Initialize arrays to store results
deltaHdiss_values = zeros(size(T));
T0_values = zeros(size(T));

for i = 1:length(T)
    % Objective function for nonlinear regression
    objectiveFunction = @(params) sum((log(Xexp) - ...
        (-params(1) / R) .* ((1 ./ T(i)) - (1 / params(2)))).^2);
    
    % Initial guess for parameters [Delta Hdiss, T0]
    initialGuess = [Hfus, mean(T)];

    % Perform nonlinear regression to find optimal parameters
    options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton', 'OptimalityTolerance', 1.0e-7);
    [optimalParams, ~] = fminunc(objectiveFunction, initialGuess, options);

    % Store results
    deltaHdiss_values(i) = optimalParams(1);
    T0_values(i) = optimalParams(2);
end

% Calculate values for the plot
deltaHdiss_over_Hfus = deltaHdiss_values / Hfus;
T0_over_Tm = T0_values / Tm;

% Plot the graph
figure;
scatter(T0_over_Tm, deltaHdiss_over_Hfus, 'filled', 'MarkerFaceColor', 'b');
xlabel('T0 / Tm');
ylabel('ΔHdiss / ΔHfus');
title('Plot of ΔHdiss/ΔHfus vs. T0/Tm');
grid on;

% Display the results
disp('Optimal Parameters at Each Temperature:');
for i = 1:length(T)
    fprintf('T = %.2f K: ΔHdiss = %.2f J/mol, T0 = %.2f K\n', T(i), deltaHdiss_values(i), T0_values(i));
end
