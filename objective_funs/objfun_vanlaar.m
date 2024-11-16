function [SSE, Xpred, ln_gamma2] = objfun_vanlaar(params, Xexp, T, Tfus, Hfus, R)
    
    ln_gamma2 = (params(1) ./ (R * T .* (1 + (params(1) .* Xexp ./ (params(2) .* (1 - Xexp))).^2)));
    
    % Calculate x_pred using SLE equation
    Xpred = exp(-(Hfus/R) * ((1./T) - (1/Tfus)) - ln_gamma2);

    % Calculate SSE
    SSE = sum((Xexp - Xpred).^2);
end
