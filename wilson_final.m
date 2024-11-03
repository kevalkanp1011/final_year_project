
clc; close all; clear all;
R = 8.314; % Gas Constant [J/mol K]


Xexp = [0.334,0.323,0.337,0.350,0.365]; % Row vector
T = 273.15 + [47.7,43.7,48.4,52.4,57.4]; % Convert to Kelvin
v2 = 96;
v1 = 85.54;

Hfus =  25.4 * 1000; % [J/mol]
Tfus =  158.7 + 273.15; % [K]

% Improved initial guess for a and b
initialGuess = [-1000,-1000];

options = optimset('Display','iter','MaxIter', 10000,'MaxFunEvals', 1.0e10,'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns',@optimplotfval);
[optimalParams] = fminsearch(@(optimalParams)objfun(optimalParams, Xexp, T ,R, Hfus, Tfus, v1, v2), initialGuess, options);

a = optimalParams(1);
b = optimalParams(2);

% Display the results
fprintf('Regressed parameter a: %.2f\n', a);
fprintf('Regressed parameter b: %.2f\n', b);

[~,Xpred, ln_gamma2] = objfun(optimalParams, Xexp, T ,R, Hfus, Tfus, v1, v2);
fprintf('X pred is: %.2f\n', Xpred);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));

T_plot = 273.15:0.1:343.15; % Temperature range in Kelvin



% Plotting
plot(T, Xpred, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');

% Customize the plot appearance
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('X(Solubility)');
title('Solubility vs. Temperature (Wilson Equation)');
grid on;
legend('Location', 'best');



function [SSE,x_pred,ln_gamma2]  = objfun(optimalParams, Xexp, T ,R, Hfus, Tfus, v1, v2)
    a = optimalParams(1);
    b = optimalParams(2);
    z1 = (v2/v1).*exp(-(a./(R.*T)));
    z2 = (v1/v2).*exp(-(b./(R.*T)));
    x1 = 1 - Xexp;
    x2 =Xexp;
    
    ln_gamma2 = -log(x2 + x1.*z2) + x1.*(z2./(x2 + x1.*z2) - z1./(x1 + x2.*z1));

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