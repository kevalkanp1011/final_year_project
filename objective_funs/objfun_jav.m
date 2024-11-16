function [SSE, Xpred] = objfun_jav(params, Xexp, T)
    A1 = params(1);
    B1 = params(2);
    A2 = params(3);
    B2 = params(4);
    J0 = params(5);
    J1 = params(6);
    J2 = params(7);

    w1 = 0.5; % Weight fraction of solvent 1 (can be adjusted)
    w2 = 1 - w1;
    
    ln_xwT = w1 * (A1 + B1 ./ T) + w2 * (A2 + B2 ./ T) + (w1 * w2 ./ T) .* (J0 + J1 * (w1 - w2) + J2 * (w1 - w2).^2);
    Xpred = exp(ln_xwT);
    
    % Calculate SSE
    SSE = sum((Xexp - Xpred).^2);
end