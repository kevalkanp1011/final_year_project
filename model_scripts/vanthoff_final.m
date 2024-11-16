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


% Initial guess for parameters [Delta Hdiss, T0]
initialGuess = [Hfus, mean(T)];

% Perform nonlinear regression to find optimal parameters
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton', 'OptimalityTolerance', 1.0e-7);

[optimalParams] = fminunc(@(params)objfun_vanthoff(params, Xexp, T, R), initialGuess, options);

% Extract the regressed parameters
deltaHdiss = optimalParams(1);
T0 = optimalParams(2);

% Calculate predicted values
[~, Xpred] = objfun_vanthoff(optimalParams, Xexp, T, R);


% Display the results
fprintf('Xpred is: %.2f\n', Xpred);
fprintf('Regressed parameter ΔHdiss: %.2f J/mol\n', deltaHdiss);
fprintf('Regressed parameter T0: %.2f K\n', T0 -273.15);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));


% Plot the results
scatter(T, Xpred,'MarkerEdgeColor', 'b', 'DisplayName', 'predicted');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'experimental');
xlim([0+273.15, 70+273.15]);
xlabel('Temperature (K)');
ylabel('X (Solubility)');
title('Solubility vs. Temperature (Vanthoff)');
grid on;
legend('Location', 'best');

deltaHdiss = optimalParams(1); 
T0 = optimalParams(2); % Calculations for the plot 
deltaHdiss_over_Hfus = deltaHdiss / Hfus; 
T0_over_Tm = T0 / Tfus;

figure; 
scatter(T0_over_Tm, deltaHdiss_over_Hfus, 'filled', 'MarkerFaceColor', 'b'); 
xlabel('T0 / Tm'); 
ylabel('ΔHdiss / ΔHfus'); 
title('Plot of ΔHdiss/ΔHfus vs. T0/Tm'); 
grid on;


