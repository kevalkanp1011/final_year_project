clc; close all; clear;

% Add path to the directory containing objfun_uniquac.m 
addpath('objective_funs');
addpath('utils');

% Load experimental data from a separate file
data = load('data.mat');
Xexp = data.Xexp;
T = data.T;
Hfus = data.Hfus;
Tfus = data.Tfus;
R = data.R;
v1 = data.v1;
v2 = data.v2;

% Initial guess for parameters [DeltaHdiss, T0] 
initialGuess = [Hfus, mean(T)];

% Optimization options
options = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval); 


% Perform the optimization
[optimalParams] = fminsearch(@(params) objfun_uniquac(params, Xexp, T, R, Hfus, Tfus, v1, v2), initialGuess, options);

% Extract optimized parameters
deltaHdiss = optimalParams(1);
T0 = optimalParams(2);


% Calculate predicted values
[~, Xpred, ln_gamma2] = objfun_uniquac(optimalParams, Xexp, T, R, Hfus, Tfus, v1, v2);


% Display the results
fprintf('Regressed parameter ΔHdiss: %.2f J/mol\n', deltaHdiss);
fprintf('Regressed parameter T0: %.2f K\n', T0);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));

% Plotting Predicted vs experimental solubility
figure;
plot(T - 273.15, Xexp, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
hold on;
plot(T_predict - 273.15, Xpred, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
xlabel('Temperature (°C)');
ylabel('X (Solubility)');
title('Solubility vs. Temperature (JAA Model)');
grid on;
legend('Location', 'best');

% Plotting ln(gamma_2) vs. Temperature
figure;
scatter(T, Xpred, 'MarkerEdgeColor', 'b', 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('ln(\gamma_2)');
title('ln(\gamma_2) vs. Temperature (UNIQUAC Model)');
grid on;
legend('Location', 'best');



