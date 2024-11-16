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

% Improved initial guess for a and b
initialGuess = [-1000, -1000];

% Optimization options
options = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval);

% Perform the optimization
[optimalParams] = fminsearch(@(optimalParams) objfun_wilson(optimalParams, Xexp, T, R, Hfus, Tfus, v1, v2), initialGuess, options);

% Extract optimized parameters
a = optimalParams(1);
b = optimalParams(2);

% Display the results
fprintf('Regressed parameter a: %.2f\n', a);
fprintf('Regressed parameter b: %.2f\n', b);

% Calculate predicted values
[~, Xpred, ln_gamma2] = objfun_wilson(optimalParams, Xexp, T, R, Hfus, Tfus, v1, v2);

% Display calculated Xpred and S
fprintf('X pred is: %.2f\n', Xpred);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));


% Plotting solubility vs temperature
figure;
scatter(T, Xpred, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility'); % Blue scatter plot
hold on;
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility'); % Red scatter plot
xlabel('Temperature (K)');
ylabel('X (Solubility)');
title('Solubility vs. Temperature (Wilson Equation)');
grid on;
legend('Location', 'best');


% Plotting gamma2 vs temperature
figure;
scatter(T, ln_gamma2, 'MarkerEdgeColor', 'b', 'DisplayName', 'Predicted ln(\gamma2)');
xlabel('Temperature (K)');
ylabel('gamma2 (-)');
title('gamma2 vs. Temperature (Wilson Equation)');
grid on;
legend('Location', 'best');
