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

% Improved initial guess for a and b
initialGuess = 0;

options = optimset('Display','iter','MaxFunEvals', 1.0e10,'MaxIter', 1e10,'TolFun', 1.0e-10, 'TolX', 1.0e-15,'PlotFcns',@optimplotfval);
[optimalParams] = fminsearch(@(params)objfun_margules(params, Xexp, T ,R, Hfus, Tfus), initialGuess, options);
[SSE,Xpred, ln_gamma2] = objfun_margules(optimalParams, Xexp, T ,R, Hfus, Tfus);

a21 = optimalParams;

% Display the results
fprintf('Regressed parameter a21: %.2f\n', a21);
fprintf('X pred is: %.2f\n', Xpred);
fprintf('SSE is: %.10f\n', SSE);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));


% Plotting
figure;
scatter(T, Xpred, 'MarkerEdgeColor', 'b', 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
xlabel('Temperature (K)');
ylabel('X(Solubility)');
title('Solubility vs. Temperature (Margules Equation)');
grid on;
legend('Location', 'best');

figure;
scatter(T, ln_gamma2, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Predicted lngamma2');
xlabel('Temperature (K)');
ylabel('gamma2(-)');
title('gamma2 vs. Temperature (Margules Equation)');
grid on;
legend('Location', 'best')





