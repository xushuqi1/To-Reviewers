function TES = smootheffscorebeforelitest(TE, ninputs, noutputs, varargin)
%   Corrected version for boundary smoothing of DEA efficiency scores (Simar & Zelenyuk, 2006)
%   Author: Adapted by User

    defopt = {0.05, 1.0}; % alpha as decimal (e.g., 0.05 for 5%)
    defopt(1:length(varargin)) = varargin;
    [alpha, effval] = defopt{:};

    k = size(TE,1);
    
    % Find DMUs with efficiency score = effval
    ind = find(abs(effval - TE) <= 1e-14);
    
    % Calculate gamma_n as per theory
    gamma_n = k^(-2/(ninputs + noutputs + 1));
    
    % Smooth efficient DMUs with positive perturbation
    TES = TE;
    TES(ind) = TE(ind) + gamma_n .* rand(length(ind),1);
end