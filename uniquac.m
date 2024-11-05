

% Constants 
R = 8.314; % Gas Constant [J/mol K] 
Hfus = 25.4 * 1000; % Heat of fusion [J/mol] 
Tfus = 158.7 + 273.15; % Fusion temperature [K] % Experimental Data Points for Solute 1 
Xexp = [0.334, 0.323, 0.337, 0.350, 0.365]; % Row vector 
T = 273.15 + [47.7, 43.7, 48.4, 52.4, 57.4]; % Convert to Kelvin % Volume ratios (assuming given data is for solutes) 
v2 = 96; 
v1 = 85.54; % Nonlinear regression options 
options = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval); 
% Initial guess for parameters [DeltaHdiss, T0] 
initialGuess = [Hfus, mean(T)];
% UNIQUAC model function 

objectiveFunction = @(params) sum((log(Xexp) - model(params, T, R, v1, v2, Hfus, Xexp)).^2);

% Perform nonlinear regression
[optimalParams] = fminsearch(objectiveFunction, initialGuess, options);

% Extract the regressed parameters
deltaHdiss = optimalParams(1);
T0 = optimalParams(2);

% Calculate the predicted solubility
Xpred = exp(model(optimalParams, T, R, v1, v2, Hfus, Xexp));

% Calculate S percent
S_percent = calculate_S(Xexp, Xpred);

% Display the results
fprintf('Regressed parameter ΔHdiss: %.2f J/mol\n', deltaHdiss);
fprintf('Regressed parameter T0: %.2f K\n', T0);
fprintf('S percent is: %.2f\n', S_percent);

% Plotting
figure;
plot(T, log(Xpred), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, log(Xexp), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('ln(\gamma_2)');
title('ln(\gamma_2) vs. Temperature (UNIQUAC Model)');
grid on;
legend('Location', 'best');

% Plotting
figure;
plot(T, Xpred, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('Xpred');
title('Xpred vs. Temperature (UNIQUAC Model)');
grid on;
legend('Location', 'best');


function lng = model(params, T, R, v1, v2,Hfus, Xexp) 
deltaHdiss = params(1); 
T0 = params(2); 
% UNIQUAC parameters 
q1 = v1 / R; 
% Segment fraction for component 1 
q2 = v2 / R;
% Segment fraction for component 2 
r1 = v1 / Hfus; % Area fraction for component 1 
r2 = v2 / Hfus; % Area fraction for component 2 
% Volume fraction 
phi1 = Xexp * r1 / (Xexp * r1 + (1 - Xexp) * r2); 
phi2 = (1 - Xexp) * r2 / (Xexp * r1 + (1 - Xexp) * r2); 
% Surface fraction 
theta1 = Xexp * q1 / (Xexp * q1 + (1 - Xexp) * q2); 
theta2 = (1 - Xexp) * q2 / (Xexp * q1 + (1 - Xexp) * q2); 
% Temperature dependent part 
tau12 = deltaHdiss ./ (R * T) .* (1 - T0 ./ T); 
tau21 = deltaHdiss ./ (R * T) .* (1 - T ./ T0); 
% Activity coefficients 
lng = q1 * log(theta1 / phi1) + q2 * log(theta2 / phi2) + phi1 * tau12 * ((1 - theta1 ./ theta2) / theta2) + phi2 * tau21 * ((1 - theta2 ./ theta1) / theta1); 
end

function S = calculate_S(x_experimental, x_calculated)
    N = length(x_experimental);
    sum_squared_errors = sum(((x_experimental - x_calculated) ./ x_experimental).^2);
    S = sqrt(1 / (N - 1) * sum_squared_errors) * 100;
end
