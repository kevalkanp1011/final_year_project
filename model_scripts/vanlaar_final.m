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


% Initial guess for parameters
initial_guess = [-1000, -1000];

% Optimization options
options = optimset('Display','iter','MaxIter', 10000,'TolFun', 1.0e-10, 'TolX', 1.0e-10,'PlotFcns',@optimplotfval);

% Perform optimization
[optimalParams] = fminsearch(@(params) objfun_vanlaar(params, Xexp, T, Tfus, Hfus, R), initial_guess, options);


[~,Xpred, ln_gamma2] = objfun_vanlaar(optimalParams, Xexp, T ,Tfus, Hfus, R);

% Extract optimal parameters
a12 = optimalParams(1);
b12 = optimalParams(2);


fprintf('Regressed parameter a12: %.2f\n', a12);
fprintf('Regressed parameter b12: %.2f\n', b12);


fprintf('X pred is: %.2f\n', Xpred);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));


% Plotting
figure;
scatter(T, Xpred, 'MarkerEdgeColor', 'b', 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
xlabel('Temperature (K)');
ylabel('X(Solubility)');
title('Solubility vs. Temperature (Vanlaar Equation)');
grid on;
legend('Location', 'best');

% Plotting
figure;
scatter(T, ln_gamma2, 'MarkerEdgeColor', 'b', 'DisplayName', 'Predicted ln(\gamma2)');
xlabel('Temperature (K)');
ylabel('gamma2(-)');
title('gamma2 vs. Temperature (Vanlaar Equation)');
grid on;
legend('Location', 'best');

