% Experimental Data Points for Solute 1
Xexp = [0.334, 0.323, 0.337, 0.350, 0.365]; % Row vector
T = 273.15 + [47.7, 43.7, 48.4, 52.4, 57.4]; % Convert to Kelvin

% Initial guess for parameters [a_A, b_A, c_A]
initialGuess = [0,0,0];

options = optimset('Display','iter','MaxFunEvals', 1.0e10,'MaxIter', 10000, 'TolFun', 1.0e-6, 'TolX', 1.0e-07,'PlotFcns',@optimplotfval);
[optimalParams] = fminsearch(@(optimalParams)objfun(optimalParams,Xexp, T), initialGuess, options);


% Extract the regressed parameters
aA = optimalParams(1);
bA = optimalParams(2);
cA = optimalParams(3);

fprintf('Regressed parameter a_A: %.2f\n', aA);
fprintf('Regressed parameter b_A: %.2f\n', bA);
fprintf('Regressed parameter c_A: %.2f\n', cA);



rp_params = [-0.5, -820, 0.35];

[F,~, Xpred] = objfun(optimalParams, Xexp, T);

fprintf('Xpred is: %.3f\n', Xpred);
fprintf('SSE is: %.10f\n', F);
fprintf('S percent is: %.6f\n', calculate_S(Xexp, Xpred));

T_celsius = 0:70;

T_Graph = T_celsius + 273.15;


Xpred_Range = exp(optimalParams(1) + optimalParams(2) ./ T_Graph + optimalParams(3) *log(T_Graph));

% Plot the results
plot(T_Graph, Xpred_Range, 'b-','LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'predicted');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'experimental');
%plot(T, Xexp, 'r-o');
xlim([0+273.15, 70+273.15]);
xlabel('Temperature (K)');
ylabel('X (Solubility)');
title('Solubility vs. Temperature (Apelblat)');
grid on;
legend('Location', 'best');


function [F, optimalParams, Xpred] = objfun(optimalParams, Xexp, T)

Xpred = exp(optimalParams(1) + optimalParams(2) ./ T + optimalParams(3) *log(T));

F = sum((Xexp - Xpred).^2);

end

function S = calculate_S(x_experimental, x_calculated)

% Calculate the number of data points
N = length(x_experimental);

% Calculate the sum of squared relative errors
sum_squared_errors = sum(((x_experimental - x_calculated) ./ x_experimental).^2);

% Calculate S
S = sqrt(1 / (N - 1) * sum_squared_errors) * 100;

end