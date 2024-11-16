function lng = uniquac_model(params, T, R, v1, v2,Hfus, Xexp) 
deltaHdiss = params(1); 
T0 = params(2); 
% UNIQUAC parameters 
q1 = v1 / R; 
% Segment fraction for component 1 
q2 = v2 / R;
% Segment fraction for component 2 
r1 = v1 / Hfus; % Area fraction for component 1 
r2 = v2 / Hfus; % Area fraction for component 2 
% Volume fraction 
phi1 = Xexp * r1 / (Xexp * r1 + (1 - Xexp) * r2); 
phi2 = (1 - Xexp) * r2 / (Xexp * r1 + (1 - Xexp) * r2); 
% Surface fraction 
theta1 = Xexp * q1 / (Xexp * q1 + (1 - Xexp) * q2); 
theta2 = (1 - Xexp) * q2 / (Xexp * q1 + (1 - Xexp) * q2); 
% Temperature dependent part 
tau12 = deltaHdiss ./ (R * T) .* (1 - T0 ./ T); 
tau21 = deltaHdiss ./ (R * T) .* (1 - T ./ T0); 
% Activity coefficients 
lng = q1 * log(theta1 / phi1) + q2 * log(theta2 / phi2) + phi1 * tau12 * ((1 - theta1 ./ theta2) / theta2) + phi2 * tau21 * ((1 - theta2 ./ theta1) / theta1); 
end