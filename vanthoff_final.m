% Experimental Data Points for Solute 1
Xexp = [0.334, 0.323, 0.337, 0.350, 0.365]; % Row vector
T = 273.15 + [47.7, 43.7, 48.4, 52.4, 57.4]; % Convert to Kelvin

Hfus =  25.4 * 1000;
% Constants
R = 8.314; % Gas constant in J/(mol*K)

% Objective function for nonlinear regression
objectiveFunction = @(params) sum((log(Xexp) - ...
    (-params(1) / R) .* ((1 ./ T) - (1 / params(2)))).^2);

% Initial guess for parameters [Delta Hdiss, T0]
initialGuess = [Hfus, mean(T)];

% Perform nonlinear regression to find optimal parameters
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton', 'OptimalityTolerance', 1.0e-7);
[optimalParams, ~] = fminunc(@(optimalParams)objfun(optimalParams, Xexp, T, R), initialGuess, options);

% Extract the regressed parameters
deltaHdiss = optimalParams(1);
T0 = optimalParams(2);

Xpred = exp((-optimalParams(1) / R) .* ((1 ./ T) - (1 / optimalParams(2))));

% Display the results
fprintf('Xpred is: %.2f\n', Xpred);
fprintf('Regressed parameter ΔHdiss: %.2f J/mol\n', deltaHdiss);
fprintf('Regressed parameter T0: %.2f K\n', T0 -273.15);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));

T_celsius = 0:70;

T_Graph = T_celsius + 273.15;


Xpred_Range = exp((-optimalParams(1) / R) .* ((1 ./ T_Graph) - (1 / optimalParams(2))));

% Plot the results
plot(T_Graph, Xpred_Range, 'b-','LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'predicted');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'experimental');
%plot(T, Xexp, 'r-o');
xlim([0+273.15, 70+273.15]);
xlabel('Temperature (K)');
ylabel('X (Solubility)');
title('Solubility vs. Temperature (Vanthoff)');
grid on;
legend('Location', 'best');


function [F, optimalParams, Xpred] = objfun(optimalParams, Xexp, T, R)

Xpred = zeros(1, length(Xexp));

for p = 1:length(Xexp)
    Xpred = exp((-optimalParams(1) / R) .* ((1 ./ T) - (1 / optimalParams(2))));
end

F = sum((log(Xexp) - ...
    (-optimalParams(1) / R) .* ((1 ./ T) - (1 / optimalParams(2)))).^2);

end

function S = calculate_S(x_experimental, x_calculated)

% Calculate the number of data points
N = length(x_experimental);

% Calculate the sum of squared relative errors
sum_squared_errors = sum(((x_experimental - x_calculated) ./ x_experimental).^2);

% Calculate S
S = sqrt(1 / (N - 1) * sum_squared_errors) * 100;

end