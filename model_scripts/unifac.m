% Experimental Data Points for Solute 1
Xexp = [0.334, 0.323, 0.337, 0.350, 0.365]; % Row vector
T = 273.15 + [47.7, 43.7, 48.4, 52.4, 57.4]; % Convert to Kelvin

% Constants
Hfus = 25.4 * 1000; % J/mol
Tfus = 158.7 + 273.15; % K
R = 8.314; % Gas constant in J/(mol*K)


% Placeholder for initial parameters (needs to be defined)
initial_params = [0.5, 0.5, 10, 1, 1, 0.5, 0.5, 1, 1, 1, 0.5, 1, 1, 0.5, 1,1];

objectiveFunction = @(params)objfun(params, Xexp,Hfus, Tfus, T, R);

options = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval); 
% Function to estimate parameters using optimization (example placeholder)
[estimated_params, ~] = fminsearch(objectiveFunction, initial_params, options);




% Display estimated parameters
disp('Estimated Parameters:');
disp(estimated_params);


function objValue = objfun(params, Xexp,Hfus, Tfus, T, R)

phi_i = params(1); theta_i = params(2); z = params(3); q_i = params(4); l_i = params(5); x_i = params(6); x_J = params(7); l_J = params(8); nu_k_i = params(9); Gamma_k = params(10); Gamma_k_i = params(11); theta_m = params(12); Psi_mk = params(13); Psi_km = params(14); theta_n = params(15); Psi_nm = params(16);

ln_gamma2 = @(phi_i, theta_i, z, q_i, l_i, x_i, x_J, l_J, nu_k_i, Gamma_k, Gamma_k_i, theta_m, Psi_mk, Psi_km, theta_n, Psi_nm) ...
    log(phi_i/x_i) + (z/2) * q_i * log(theta_i/phi_i) + l_i - (phi_i/x_i) * sum(x_J * l_J) + ...
    sum(nu_k_i * (log(Gamma_k) - log(Gamma_k_i)));
ln_gamma2_val = ln_gamma2(phi_i, theta_i, z, q_i, l_i, x_i, x_J, l_J, nu_k_i, Gamma_k, Gamma_k_i, theta_m, Psi_mk, Psi_km, theta_n, Psi_nm);

Xpred = @(T, ln_gamma2) exp(-(Hfus/R) * ((1./T) - (1/Tfus)) - ln_gamma2);

objValue = sum((Xexp - Xpred(T, ln_gamma2_val)).^2);

end
