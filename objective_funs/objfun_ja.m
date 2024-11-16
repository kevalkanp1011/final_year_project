function [SSE, Xpred] = objfun_ja(optimalParams, Xexp, T)
    J0 = optimalParams(1);
    J1 = optimalParams(2);
    J2 = optimalParams(3);
    
    w1 = 0.5; % Example weight fraction, adjust as needed
    w2 = 1 - w1;
    
    ln_xwT = w1 * log(Xexp) + w2 * log(Xexp) + (w1 * w2 ./ T) .* (J0 + J1 * (w1 - w2) + J2 * (w1 - w2).^2);
    Xpred = exp(ln_xwT);

    % Calculate SSE
    SSE = sum((Xexp - Xpred).^2);
end