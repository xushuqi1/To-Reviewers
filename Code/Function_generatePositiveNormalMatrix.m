function P = Function_generatePositiveNormalMatrix(n, n1, n2)    %   generate data with normal distribution
    % Initialize the output matrix P
    P = zeros(n, n1 + n2);
    
    % Parameters for the normal distribution
    mu = 1000;  % Mean
    sigma = sqrt(100);  % Standard deviation, since the variance is 100

    % Fill the matrix with values
    for col = 1:(n1 + n2)
        % Repeat until all values in the column are greater than 0
        while true
            % Generate a column vector of n values from N(1000, 100)
            colValues = normrnd(mu, sigma, [n, 1]);
            
            % Check if all values in the column are greater than 0
            if all(colValues > 0)
                P(:, col) = colValues;  % Assign the column vector to the matrix
                break;  % Exit the loop if all values are greater than 0
            end
        end
    end
end
