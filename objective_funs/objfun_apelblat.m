function [SSE, Xpred] = objfun_apelblat(params, Xexp, T)

Xpred = exp(params(1) + params(2) ./ T + params(3) *log(T));

SSE = sum((Xexp - Xpred).^2);

end