clc; close all; clear;

% Load experimental data from a separate file
data = load('data.mat');
Xexp = data.Xexp;
T = data.T;
Hfus = data.Hfus;
Tfus = data.Tfus;
R = data.R;
v1 = data.v1;
v2 = data.v2;

% Van't Hoff Model
initialGuess_vanthoff = [Hfus, mean(T)];
options_vanthoff = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton', 'OptimalityTolerance', 1.0e-7);
[optimalParams_vanthoff] = fminunc(@(params)objfun_vanthoff(params, Xexp, T, R), initialGuess_vanthoff, options_vanthoff);
[~, Xpred_vanthoff] = objfun_vanthoff(optimalParams_vanthoff, Xexp, T, R);

% LamdaH Model
initialGuess_lamdah = [0,0];
options_lamdah = optimset('Display','iter','TolFun', 1.0e-10, 'TolX', 1.0e-05,'PlotFcns',@optimplotfval);
[optimalParams_lamdah] = fminsearch(@(params)objfun_lamdah(params,Xexp, T, Tfus), initialGuess_lamdah, options_lamdah);
[~, Xpred_lamdah] = objfun_lamdah(optimalParams_lamdah,Xexp, T, Tfus);

% Apelblat Model
initialGuess_apelblat = [0,0,0];
options_apelbalt = optimset('Display','iter','MaxFunEvals', 1.0e10,'MaxIter', 10000, 'TolFun', 1.0e-6, 'TolX', 1.0e-07,'PlotFcns',@optimplotfval);
[optimalParams_apelblat] = fminsearch(@(optimalParams)objfun_apelblat(optimalParams,Xexp, T), initialGuess_apelblat, options_apelbalt);
[~,Xpred_apelblat] = objfun_apelblat(optimalParams_apelblat, Xexp, T);

% JA Model
initialGuess_ja = [-1000, -1000, -1000];
options_ja = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval);
[optimalParams_ja] = fminsearch(@(optimalParams) objfun_ja(optimalParams, Xexp, T), initialGuess_ja, options_ja);
[~, Xpred_ja] = objfun_ja(optimalParams_ja, Xexp, T);

% JAV Model
initialGuess_jav = [0, 0, 0, 0, 0, 0, 0];
options_jav = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval);
[optimalParams_jav] = fminsearch(@(params) objfun_jav(params, Xexp, T), initialGuess_jav, options_jav);
[~, Xpred_jav] = objfun_jav(optimalParams_jav, Xexp, T);

% JAA Model
initialGuess_jaa = [1, 1, 1, 1, 1, 1, 1, 1, 1];
options_jaa = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval);
[optimalParams_jaa] = fminsearch(@(params) objfun_jaa(params, Xexp, T), initialGuess_jaa, options_jaa);
[~, Xpred_jaa] = objfun_jaa(optimalParams_jaa, Xexp, T);



% Margules model
initial_guess_margules = 0;
options_margules = optimset('Display','iter','MaxIter', 10000,'TolFun', 1.0e-10, 'TolX', 1.0e-10,'PlotFcns',@optimplotfval);
[optimalParams_margules] = fminsearch(@(params)objfun_margules(params, Xexp, T ,R, Hfus, Tfus), initial_guess_margules, options_margules);
[~,Xpred_margules, ln_gamma2_margules] = objfun_margules(optimalParams_margules, Xexp, T ,R, Hfus, Tfus);

% Vanlaar model
initial_guess_vanlaar = [-1000, -1000];
options_vanlaar = optimset('Display','iter','MaxIter', 10000,'TolFun', 1.0e-10, 'TolX', 1.0e-10,'PlotFcns',@optimplotfval);
[optimalParams_vanlaar] = fminsearch(@(params) objfun_vanlaar(params, Xexp, T, Tfus, Hfus, R), initial_guess_vanlaar, options_vanlaar);
[~,Xpred_vanlaar, ln_gamma2_vanlaar] = objfun_vanlaar(optimalParams_vanlaar, Xexp, T ,Tfus, Hfus, R);

% Wilson Model
initialGuess_wilson = [-1000, -1000];
options_wilson = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval);
[optimalParams_wilson] = fminsearch(@(params) objfun_wilson(params, Xexp, T, R, Hfus, Tfus, v1, v2), initialGuess_wilson, options_wilson);
[~, Xpred_wilson, ln_gamma2_wilson] = objfun_wilson(optimalParams_wilson, Xexp, T, R, Hfus, Tfus, v1, v2);

% NRTL Model
initialGuess_nrtl = [-1000,-1000];
options_nrtl = optimset('Display','iter','TolFun', 1.0e-10, 'TolX', 1.0e-05,'PlotFcns',@optimplotfval);
[optimalParams_nrtl] = fminsearch(@(params)objfun_nrtl(params, Xexp, T, Hfus, Tfus, R), initialGuess_nrtl, options_nrtl);
[~,Xpred_nrtl, ln_gamma2_nrtl] = objfun_nrtl(optimalParams_nrtl, Xexp, T, Hfus, Tfus, R);


% UNIQUAC Model
initialGuess_uniquac = [Hfus, mean(T)];
options_uniquac = optimset('Display', 'iter', 'MaxIter', 10000, 'MaxFunEvals', 1.0e10, 'TolFun', 1.0e-10, 'TolX', 1.0e-07, 'PlotFcns', @optimplotfval);
[optimalParams_uniquac] = fminsearch(@(params) objfun_uniquac(params, Xexp, T, R, Hfus, Tfus, v1, v2), initialGuess_uniquac, options_uniquac);
[~, Xpred_uniquac, ln_gamma2_uniquac] = objfun_uniquac(optimalParams_uniquac, Xexp, T, R, Hfus, Tfus, v1, v2);

% Plotting combined results
figure;
hold on;
scatter(T, Xexp, 'filled', 'MarkerFaceColor', 'r','DisplayName', 'Experimental Data'); % Red scatter plot
scatter(T, Xpred_vanthoff, 'MarkerEdgeColor', 'b','LineWidth', 1.5,'DisplayName', 'Vant Hoff'); % Blue scatter plot
scatter(T, Xpred_apelblat, 'MarkerEdgeColor', 'c', 'LineWidth', 1.5,'DisplayName', 'Apelblat'); % Cyan scatter plot
scatter(T, Xpred_wilson, 'MarkerEdgeColor', 'g', 'LineWidth', 1.5,'DisplayName', 'Wilson'); % Green scatter plot
scatter(T, Xpred_uniquac, 'MarkerEdgeColor', 'm', 'LineWidth', 1.5,'DisplayName', 'UNIQUAC'); % Magenta scatter plot
scatter(T, Xpred_vanlaar, 'MarkerEdgeColor', 'y', 'LineWidth', 1.5,'DisplayName', 'Vanlaar'); % Yellow scatter plot
scatter(T, Xpred_nrtl, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5,'DisplayName', 'NRTL'); % Black scatter plot
scatter(T, Xpred_margules, 'MarkerEdgeColor', '#800080', 'LineWidth', 1.5,'DisplayName', 'Margules'); % Purple scatter plot
scatter(T, Xpred_lamdah, 'MarkerEdgeColor', '#FFA500', 'LineWidth', 1.5,'DisplayName', 'Lamda H'); % Purple scatter plot
scatter(T, Xpred_ja, 'MarkerEdgeColor', '#008080', 'LineWidth', 1.5,'DisplayName', 'Jouyban Acree'); % Purple scatter plot
scatter(T, Xpred_jav, 'MarkerEdgeColor', '#FFC0CB', 'LineWidth', 1.5,'DisplayName', 'JAV'); % Purple scatter plot
scatter(T, Xpred_jaa, 'MarkerEdgeColor', '#A52A2A', 'LineWidth', 1.5,'DisplayName', 'JAA'); % Purple scatter plot
xlabel('Temperature (K)');
ylabel('X (Solubility)');
title('Solubility vs. Temperature (Combined Models)');
grid on;
legend('Location', 'best');

% Plotting ln(gamma_2) 
figure;
hold on;
scatter(T, ln_gamma2_margules,'filled', 'MarkerEdgeColor', 'g', 'DisplayName', 'Margules');
scatter(T, ln_gamma2_vanlaar,'MarkerEdgeColor', 'b', 'DisplayName', 'Vanlaar');
scatter(T, ln_gamma2_wilson, 'MarkerEdgeColor', 'g', 'DisplayName', 'Wilson');
scatter(T, ln_gamma2_nrtl, 'MarkerEdgeColor', 'c', 'DisplayName', 'NRTL');
scatter(T, ln_gamma2_uniquac, 'MarkerEdgeColor', 'm', 'DisplayName', 'UNIQUAC');
xlabel('Temperature (K)');
ylabel('ln(\gamma_2)');
title('ln(\gamma_2) vs. Temperature (Vanlaar, Wilson, NRTL and UNIQUAC Models)');
grid on;
legend('Location', 'best');
