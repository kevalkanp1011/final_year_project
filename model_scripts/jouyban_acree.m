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

% Initial guess for Jouyban-Acree parameters (J0, J1, J2)
initialGuess = [-1000, -1000, -1000];

options = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval);
[optimalParams] = fminsearch(@(optimalParams) objfun_ja(optimalParams, Xexp, T), initialGuess, options);

[~, Xpred] = objfun_ja(optimalParams, Xexp, T);

J0 = optimalParams(1);
J1 = optimalParams(2);
J2 = optimalParams(3);

% Display the results
fprintf('Regressed parameter J0: %.2f\n', J0);
fprintf('Regressed parameter J1: %.2f\n', J1);
fprintf('Regressed parameter J2: %.2f\n', J2);


fprintf('X pred is: %.3f\n', Xpred);
fprintf('S percent is: %.7f\n', calculate_S(Xexp, Xpred));

% Plotting predicted vs experimental solubility
figure;
scatter(T, Xpred, 'MarkerFaceColor', 'b','LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
hold on;
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
xlim([min(T), max(T)]);
xlabel('Temperature (K)');
ylabel('X (Solubility)');
title('Solubility vs. Temperature (Jouyban-Acree Model)');
grid on;
legend('Location', 'best');



