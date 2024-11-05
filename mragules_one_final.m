% Initialization
clc; close all; clear all;
R = 8.314; % Gas Constant [J/mol K]


% Experimental Data Points for Solute 1
Xexp = [0.334,0.323,0.337,0.350,0.365]; % Row vector
T = 273.15 + [47.7,43.7,48.4,52.4,57.4]; % Convert to Kelvin

% Results of [g12, g21 in KJ/mol] is [ 1.18 -5.02]

% % 
% Constants for Solute 1
Hfus =  25.4 * 1000; % [J/mol]
Tfus =  158.7 + 273.15; % [K]

% Improved initial guess for a and b
initialGuess = 0;

options = optimset('Display','iter','MaxFunEvals', 1.0e10,'MaxIter', 1e10,'TolFun', 1.0e-10, 'TolX', 1.0e-15,'PlotFcns',@optimplotfval);
[optimalParams] = fminsearch(@(optimalParams)objfun(optimalParams, Xexp, T ,R, Hfus, Tfus), initialGuess, options);

a21 = optimalParams;
% Display the results
fprintf('Regressed parameter a21: %.2f\n', a21);

rp_params = -8030;
[SSE,Xpred, ln_gamma2] = objfun(optimalParams, Xexp, T ,R, Hfus, Tfus);
fprintf('X pred is: %.2f\n', Xpred);
fprintf('SSE is: %.10f\n', SSE);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));

T_plot = 273.15:0.1:343.15; % Temperature range in Kelvin



% Plotting
figure;
plot(T, Xpred, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');

% Customize the plot appearance
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('X(Solubility)');
title('Solubility vs. Temperature (Margules Equation)');
grid on;
legend('Location', 'best');

figure;
plot(T, ln_gamma2, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, ln_gamma2, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
% Customize the plot appearance
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('gamma2(-)');
title('gamma2 vs. Temperature (Margules Equation)');
grid on;
legend('Location', 'best')




function [SSE,x_pred, ln_gamma2]  = objfun(optimalParams, Xexp, T ,R, Hfus, Tfus)
    A21 = optimalParams./(R*T);

    % Calculate predicted activity coefficients using Margules equation
    ln_gamma2 = A21 .* (1-Xexp).^2;

    % Calculate x_pred using SLE equation
    x_pred = exp(-(Hfus/R) * ((1./T) - (1/Tfus)) - ln_gamma2);

    % Calculate SSE
    SSE = sum((Xexp - x_pred).^2);
end

function S = calculate_S(x_experimental, x_calculated)

% Calculate the number of data points
N = length(x_experimental);

% Calculate the sum of squared relative errors
sum_squared_errors = sum(((x_experimental - x_calculated) ./ x_experimental).^2);

% Calculate S
S = sqrt(1 / (N - 1) * sum_squared_errors) * 100;

end