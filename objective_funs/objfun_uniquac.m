function [SSE, Xpred, ln_gamma2] = objfun_uniquac(params, Xexp, T, R, Hfus, Tfus, v1, v2)
   

    addpath('./model_funs');
    ln_gamma2 = uniquac_model(params, T, R, v1, v2, Hfus, Xexp);

    % Calculate x_pred using SLE equation
    Xpred = exp(-(Hfus/R) * ((1 ./ T) - (1/Tfus)) - ln_gamma2);

    % Calculate SSE
    SSE = sum((Xexp - Xpred).^2);
end