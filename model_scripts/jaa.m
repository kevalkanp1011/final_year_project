clc; close all; clear;

% Add path to the directory containing objfun_uniquac.m 
addpath('objective_funs');
addpath('utils');

% Load experimental data from a separate file
data = load('data.mat');
Xexp = data.Xexp;
T = data.T;

% Initial guess for model parameters [A1, B1, C1, A2, B2, C2, J0, J1, J2]
initialGuess = [1, 1, 1, 1, 1, 1, 1, 1, 1];

% Perform parameter estimation using fminsearch
options = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval);
[optimalParams] = fminsearch(@(params) objfun_jaa(params, Xexp, T), initialGuess, options);
[~, Xpred] = objfun_jaa(optimalParams, Xexp, T);

% Extract estimated parameters
A1 = optimalParams(1);
B1 = optimalParams(2);
C1 = optimalParams(3);
A2 = optimalParams(4);
B2 = optimalParams(5);
C2 = optimalParams(6);
J0 = optimalParams(7);
J1 = optimalParams(8);
J2 = optimalParams(9);

% Display the results
fprintf('Estimated parameters for JAA model:%.4f\n', [A1, B1, C1, A2, B2, C2, J0, J1, J2]);
fprintf('S percent is: %.7f\n', calculate_S(Xexp, Xpred));


% Plotting predicted vs experimental solubility
figure;
plot(T - 273.15, Xexp, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
hold on;
plot(T - 273.15, Xpred, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
xlabel('Temperature (°C)');
ylabel('X (Solubility)');
title('Solubility vs. Temperature (JAA Model)');
grid on;
legend('Location', 'best');

