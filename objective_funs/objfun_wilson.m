function [SSE, Xpred, ln_gamma2] = objfun_wilson(optimalParams, Xexp, T, R, Hfus, Tfus, v1, v2)
    a = optimalParams(1);
    b = optimalParams(2);
    z1 = (v2/v1) .* exp(-(a ./ (R .* T)));
    z2 = (v1/v2) .* exp(-(b ./ (R .* T)));
    x1 = 1 - Xexp;
    x2 = Xexp;

    ln_gamma2 = -log(x2 + x1 .* z2) + x1 .* (z2 ./ (x2 + x1 .* z2) - z1 ./ (x1 + x2 .* z1));

    % Calculate x_pred using SLE equation
    Xpred = exp(-(Hfus/R) * ((1 ./ T) - (1/Tfus)) - ln_gamma2);

    % Calculate SSE
    SSE = sum((Xexp - Xpred).^2);
end
