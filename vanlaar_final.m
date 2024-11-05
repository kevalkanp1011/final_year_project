% Experimental Data Points for Solute 1
Xexp = [0.334, 0.323, 0.337, 0.350, 0.365]; % Row vector
T = 273.15 + [47.7, 43.7, 48.4, 52.4, 57.4]; % Convert to Kelvin

% Constants for Solute 1
Hfus =  25.4 * 1000; % [J/mol]
Tfus =  158.7 + 273.15; % [K]
R = 8.314; % Gas constant in J/(mol*K)

% Initial guess for parameters
initial_guess = [-1000, -1000];

% Optimization options
options = optimset('Display','iter','MaxIter', 10000,'TolFun', 1.0e-10, 'TolX', 1.0e-10,'PlotFcns',@optimplotfval);

% Perform optimization
[optimalParams] = fminsearch(@(optimalParams) objfun(optimalParams, Xexp, T, Tfus, Hfus, R), initial_guess, options);


rp_params = [-7.27*1000, -11.7*1000];
[SSE,Xpred] = objfun(optimalParams, Xexp, T ,Tfus, Hfus, R);
% Extract optimal parameters
a12 = optimalParams(1);
b12 = optimalParams(2);


fprintf('Regressed parameter a12: %.2f\n', a12);
fprintf('Regressed parameter b12: %.2f\n', b12);


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
title('Solubility vs. Temperature (Vanlaar Equation)');
grid on;
legend('Location', 'best');

% Plotting
figure;
plot(T, ln_gamma2, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Predicted Solubility');
hold on
scatter(T, ln_gamma2, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Experimental Solubility');
% Customize the plot appearance
xlim([273.15, 343.15]);
xlabel('Temperature (K)');
ylabel('gamma2(-)');
title('gamma2 vs. Temperature (Vanlaar Equation)');
grid on;
legend('Location', 'best');


function [SSE, x_pred] = objfun(optimalParams, Xexp, T, Tfus, Hfus, R)
    
    ln_gamma2 = (optimalParams(1) ./ (R * T .* (1 + (optimalParams(1) .* Xexp ./ (optimalParams(2) .* (1 - Xexp))).^2)));
    
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