function [SSE, Xpred] = objfun_vanthoff(params, Xexp, T, R)

Xpred = exp((-params(1) / R) .* ((1 ./ T) - (1 / params(2))));

% Calculate SSE
SSE = sum((Xexp - Xpred).^2);

end