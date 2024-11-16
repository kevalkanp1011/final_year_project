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

% Initial guess for parameters [lambda, h]
initialGuess = [0,0];


options = optimset('Display','iter','TolFun', 1.0e-10, 'TolX', 1.0e-05,'PlotFcns',@optimplotfval);
[optimalParams] = fminsearch(@(optimalParams)objfun_lamdah(optimalParams,Xexp, T, Tfus), initialGuess, options);

[SSE, Xpred] = objfun_lamdah(optimalParams,Xexp, T, Tfus);

lambda = optimalParams(1);
h = optimalParams(2);

% Display the results
fprintf('Regressed parameter lamda: %.3f\n', lambda);
fprintf('Regressed parameter h: %.2f\n', h);

fprintf('SSE is: %.10f\n', SSE);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));


% Plot the results
scatter(T, Xexp,  'MarkerEdgeColor', 'b','DisplayName', 'predicted');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'experimental');
xlim([0+273.15, 70+273.15]);
xlabel('Temperature (K)');
ylabel('X(Solubility)');
title('Solubility vs. Temperature (LamdaH)');
grid on;
legend('Location', 'best');




