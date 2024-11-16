function [SSE, Xpred, ln_gamma2] = objfun_nrtl(optimalParams, Xexp, T, Hfus, Tfus, R)
    % Extract parameters from the input vector
    tau12 = optimalParams(1)./(R.*T);
    tau21 = optimalParams(2)./(R.*T);

    % Calculate the non-randomness parameters (assuming alpha = 0.3)
    G12 = exp(-0.3*tau12);
    G21 = exp(-0.3*tau21);

    x1 = 1 - Xexp;
    x2 =Xexp;

    ln_gamma2 = (x1.^2) .* ((tau12.*(G12./(x2+x1.*G12)).^2) + (tau21 .* (G21 ./ ( (x1 + G21.*x2).^2) ) ));

    Xpred = exp(-(Hfus/R) * ((1./T) - (1/Tfus)) - ln_gamma2);
    
    SSE = sum((Xexp - Xpred).^2);

end