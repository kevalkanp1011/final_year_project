function [SSE, Xpred] = objfun_jaa(params, Xexp, T)
    A1 = params(1);
    B1 = params(2);
    C1 = params(3);
    A2 = params(4);
    B2 = params(5);
    C2 = params(6);
    J0 = params(7);
    J1 = params(8);
    J2 = params(9);

    w1 = 0.5; % Example weight fraction, can be adjusted
    w2 = 1 - w1;
    
    ln_xwT = w1 * (A1 + B1 ./ T + C1 .* log(T)) + w2 * (A2 + B2 ./ T + C2 .* log(T)) + (w1 * w2 ./ T) .* (J0 + J1 * (w1 - w2) + J2 * (w1 - w2).^2);
    Xpred = exp(ln_xwT);
    
    % Calculate SSE
    SSE = sum((Xexp - Xpred).^2);
end