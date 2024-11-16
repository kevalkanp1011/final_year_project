function [SSE, Xpred] = objfun_lamdah(params, Xexp, T, Tm)
    lambda = params(1);
    h = (2);

    
    Xpred = lambda ./ (exp(lambda*h*((1./T) - (1/Tm))) + lambda);


    SSE = sum((Xexp - Xpred).^2);
end
