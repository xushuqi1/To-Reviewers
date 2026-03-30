function [output] = data_generating_one_pollutinginput(n, c, s, scale_parameter) 
% Input: sample size n, c is the correlation coefficient between outputs, 
% s is the variance of the noise distribution, here we set it to zero. 
% The last column of the output is the true efficiency score.

X = 1 + 3 * rand(n, 4);                     % Generate an input matrix with a uniform distribution on [1,4]
Y = zeros(n, 1);
V = zeros(n, 1);
score = zeros(n, 1);                           

mu_1 = [0 0];                               % Simulation of inefficiency coefficients
Sigma_1 = [0.25, 0.25*c; 0.25*c, 0.25];      % Covariance matrix (check before running)
samples_1 = mvnrnd(mu_1, Sigma_1, n);        % Generate n samples
samples_1 = abs(samples_1);                  % Ensure positive values
inefficiency_desirable = samples_1(:, 1);
inefficiency_undesirable = samples_1(:, 2);

mu_2 = [0 0];                               % Simulation of noise coefficients
Sigma_2 = [s, 0; 0, s];                      % Covariance matrix (check before running)
samples_2 = mvnrnd(mu_2, Sigma_2, n);        % Generate n samples
noise_desirable = samples_2(:, 1);
noise_undesirable = samples_2(:, 2);

% Generate desirable output vector
for k=1:n
    a = X(k, :);
    Y(k) = exp(scale_parameter*log(a(1)) + scale_parameter*log(a(2)) + scale_parameter*log(a(3)) + scale_parameter*log(a(4)) + ...
           0.5*[0.1*log(a(1))*log(a(1)) - 0.1*log(a(2))*log(a(1)) - 0.1*log(a(1))*log(a(2)) + 0.1*log(a(2))*log(a(2)) + ...
           0.1*log(a(3))*log(a(3)) - 0.1*log(a(4))*log(a(3)) - 0.1*log(a(3))*log(a(4)) + 0.1*log(a(4))*log(a(4))] - ...
           inefficiency_desirable(k) + noise_desirable(k));
end                                         

% Generate undesirable output vector
for k=1:n
    a = X(k, :);
    V(k) = exp(4*scale_parameter*log(a(1)) + 0.5*0.1*log(a(1))*log(a(1)) + ...
           inefficiency_undesirable(k) + noise_undesirable(k));
end                                         

% Generate efficiency scores
for k=1:n
    score(k) = exp(-(inefficiency_desirable(k) + inefficiency_undesirable(k)) / 2);
end                                         

[output] = [X, Y, V, score];

end
