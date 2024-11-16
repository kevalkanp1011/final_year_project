function [SSE,Xpred, ln_gamma2]  = objfun_margules(params, Xexp, T ,R, Hfus, Tfus)
    A21 = params./(R*T);

    % Calculate predicted activity coefficients using Margules equation
    ln_gamma2 = A21 .* (1-Xexp).^2;

    % Calculate x_pred using SLE equation
    Xpred = exp(-(Hfus/R) * ((1./T) - (1/Tfus)) - ln_gamma2);

    % Calculate SSE
    SSE = sum((Xexp - Xpred).^2);
end