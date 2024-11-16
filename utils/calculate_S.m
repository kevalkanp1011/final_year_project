function S = calculate_S(x_experimental, x_calculated)
    % Calculate the number of data points
    N = length(x_experimental);

    % Calculate the sum of squared relative errors
    sum_squared_errors = sum(((x_experimental - x_calculated) ./ x_experimental).^2);

    % Calculate S
    S = sqrt(1 / (N - 1) * sum_squared_errors) * 100;
end
