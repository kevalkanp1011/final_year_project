% Experimental Data Points for Solute 1
Xexp = [0.334, 0.323, 0.337, 0.350, 0.365]; % Row vector
T = 273.15 + [47.7, 43.7, 48.4, 52.4, 57.4]; % Convert to Kelvin

% Constants for Solute 1
Hfus =  25.4 * 1000; % [J/mol]
Tfus =  158.7 + 273.15; % [K]
R = 8.314; % Gas constant in J/(mol*K)

% Initial guess for parameters [lambda, h]
initialGuess = [0,0];


options = optimset('Display','iter','TolFun', 1.0e-10, 'TolX', 1.0e-05,'PlotFcns',@optimplotfval);
[optimalParams] = fminsearch(@(optimalParams)objfun(optimalParams,Xexp, T, Tfus), initialGuess, options);


lambda = optimalParams(1);
h = optimalParams(2);
% Display the results
fprintf('Regressed parameter lamda: %.3f\n', lambda);
fprintf('Regressed parameter h: %.2f\n', h);

rp_params = [0.11, 2240];

[SSE, Xpred] = objfun(optimalParams,Xexp, T, Tfus);
fprintf('SSE is: %.10f\n', SSE);
fprintf('S percent is: %.2f\n', calculate_S(Xexp, Xpred));

T_celsius = 0:70;

T_Graph = T_celsius + 273.15;


Xpred_Range = lambda ./ (exp(lambda*h*((1./T_Graph) - (1/Tfus))) + lambda);

% Plot the results
plot(T_Graph, Xpred_Range, 'b-','LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'predicted');
hold on
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'experimental');
%plot(T, Xexp, 'r-o');
xlim([0+273.15, 70+273.15]);
xlabel('Temperature (K)');
ylabel('X(Solubility)');
title('Solubility vs. Temperature (LamdaH)');
grid on;
legend('Location', 'best');



function [SSE, Xpred] = objfun(optimalParams, Xexp, T, Tm)
    lambda = optimalParams(1);
    h = optimalParams(2);

    
Xpred = lambda ./ (exp(lambda*h*((1./T) - (1/Tm))) + lambda);


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