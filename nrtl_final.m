clc; close all; clear all;
R = 8.314; % Gas Constant [J/mol K]


% Experimental Data Points for Solute 1
Xexp = [0.334,0.323,0.337,0.350,0.365]; % Row vector
T = 273.15 + [47.7,43.7,48.4,52.4,57.4]; % Convert to Kelvin
Hfus =  25.4 * 1000; % [J/mol]
Tfus =  158.7 + 273.15; % [K]

initialGuess = [-1000,-1000];

options = optimset('Display','iter','TolFun', 1.0e-10, 'TolX', 1.0e-05,'PlotFcns',@optimplotfval);
[optimalParams] = fminsearch(@(optimalParams)objfun(optimalParams, Xexp, T, Hfus, Tfus, R), initialGuess, options);

a = optimalParams(1);
b = optimalParams(2);

% Display the results
fprintf('Regressed parameter a: %.2f\n', a);
fprintf('Regressed parameter b: %.2f\n', b);

[~,Xpred, ln_gamma2] = objfun(optimalParams, Xexp, T, Hfus, Tfus, R);
fprintf('X pred is: %.2f\n', Xpred);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));

T_plot = 273.15:0.1:343.15; % Temperature range in Kelvin



% Plotting
figure;
plot(T, ln_gamma2, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, ln_gamma2, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
% Customize the plot appearance
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('gamma2(-)');
title('gamma2 vs. Temperature (NRTL Equation)');
grid on;
legend('Location', 'best');


% Plotting
figure;
plot(T, Xpred, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');

% Customize the plot appearance
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('X(Solubility)');
title('Solubility vs. Temperature (NRTL Equation)');
grid on;
legend('Location', 'best');


% Objective function
function [SSE, Xpred, ln_gamma2] = objfun(optimalParams, Xexp, T, Hfus, Tfus, R)
    % Extract parameters from the input vector
    tau12 = optimalParams(1)./(R.*T);
    tau21 = optimalParams(2)./(R.*T);

    % Calculate the non-randomness parameters (assuming alpha = 0.3)
    G12 = exp(-0.3*tau12);
    G21 = exp(-0.3*tau21);

    % Calculate predicted activity coefficients using the NRTL model
    %gamma_pred = NRTL_binary(xexp, tau12, tau21, G12, G21);

    x1 = 1 - Xexp;
    x2 =Xexp;

    ln_gamma2 = (x1.^2) .* ((tau12.*(G12./(x2+x1.*G12)).^2) + (tau21 .* (G21 ./ ( (x1 + G21.*x2).^2) ) ));

    Xpred = exp(-(Hfus/R) * ((1./T) - (1/Tfus)) - ln_gamma2);
    
    SSE = sum((Xexp - Xpred).^2);

end




function S = calculate_S(x_experimental, x_calculated)

% Calculate the number of data points
N = length(x_experimental);

% Calculate the sum of squared relative errors
sum_squared_errors = sum(((x_experimental - x_calculated) ./ x_experimental).^2);

% Calculate S
S = sqrt(1 / (N - 1) * sum_squared_errors) * 100;

end