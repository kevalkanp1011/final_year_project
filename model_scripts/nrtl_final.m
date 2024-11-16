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

initialGuess = [-1000,-1000];

options = optimset('Display','iter','TolFun', 1.0e-10, 'TolX', 1.0e-05,'PlotFcns',@optimplotfval);
[optimalParams] = fminsearch(@(optimalParams)objfun_nrtl(optimalParams, Xexp, T, Hfus, Tfus, R), initialGuess, options);

[~,Xpred, ln_gamma2] = objfun_nrtl(optimalParams, Xexp, T, Hfus, Tfus, R);

a = optimalParams(1);
b = optimalParams(2);

% Display the results
fprintf('Regressed parameter a: %.2f\n', a);
fprintf('Regressed parameter b: %.2f\n', b);

fprintf('X pred is: %.2f\n', Xpred);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));

% Plotting
figure;
scatter(T, Xpred, 'MarkerEdgeColor', 'b', 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('X(Solubility)');
title('Solubility vs. Temperature (NRTL Equation)');
grid on;
legend('Location', 'best');

% Plotting
figure;
scatter(T, ln_gamma2, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'ln_gamma2');
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('ln_gamma2(-)');
title('ln_gamma2 vs. Temperature (NRTL Equation)');
grid on;
legend('Location', 'best');
