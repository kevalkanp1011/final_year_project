clc; close all; clear;

% Add path to the directory containing objfun_uniquac.m 
addpath('objective_funs');
addpath('utils');

% Load experimental data from a separate file
data = load('data.mat');
Xexp = data.Xexp;
T = data.T;
R = data.R;


% Initial guess for parameters [a_A, b_A, c_A]
initialGuess = [0,0,0];

options = optimset('Display','iter','MaxFunEvals', 1.0e10,'MaxIter', 10000, 'TolFun', 1.0e-6, 'TolX', 1.0e-07,'PlotFcns',@optimplotfval);
[optimalParams] = fminsearch(@(optimalParams)objfun_apelblat(optimalParams,Xexp, T), initialGuess, options);


% Extract the regressed parameters
aA = optimalParams(1);
bA = optimalParams(2);
cA = optimalParams(3);

fprintf('Regressed parameter a_A: %.2f\n', aA);
fprintf('Regressed parameter b_A: %.2f\n', bA);
fprintf('Regressed parameter c_A: %.2f\n', cA);


[~,Xpred] = objfun_apelblat(optimalParams, Xexp, T);

fprintf('Xpred is: %.3f\n', Xpred);
fprintf('S percent is: %.6f\n', calculate_S(Xexp, Xpred));


% Plot the results
scatter(T, Xpred,'MarkerEdgeColor', 'b',  'DisplayName', 'predicted');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'experimental');
xlim([0+273.15, 70+273.15]);
xlabel('Temperature (K)');
ylabel('X (Solubility)');
title('Solubility vs. Temperature (Apelblat)');
grid on;
legend('Location', 'best');



