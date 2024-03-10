function [output] = First_phase_twoD_Parallel(P, xx, yy) %  P refers to a n*(xx+yy) matrix.
% The input-oriented model is used here.
% To avoid the effect of the precision setting, the linear programming we use in determining the interior
% points comes from Dula (2011)
    twoD_subsample = twoD(P, xx, yy);
    n = size(twoD_subsample, 1);    
    nnn = size(P, 1);
    X1 = twoD_subsample(:, 1:xx);
    Y1 = twoD_subsample(:, xx+1:xx+yy);

    unsolvableOrNegativeK = [];
    optimalValues = zeros(1, nnn);

    parfor k = 1:nnn
        c = [zeros(n,1); 1];
        A = [X1', -ones(xx,1); -Y1', -ones(yy,1)];
        b = [P(k,1:xx)'; -P(k,xx+1:xx+yy)'];
        Aeq = [ones(1,n), 0];
        beq = 1;
        lb = [zeros(n+1,1)];
        ub = [];

        try
            [sol, fval] = linprog(c, A, b, Aeq, beq, lb, ub);
            if fval > 1e-32
                unsolvableOrNegativeK = [unsolvableOrNegativeK, k];
            end
            optimalValues(k) = fval;
        catch
            unsolvableOrNegativeK = [unsolvableOrNegativeK, k];
        end
    end

    if ~isempty(unsolvableOrNegativeK)
        unsolvableOrNegativeK = unsolvableOrNegativeK';
        exterior = P(unsolvableOrNegativeK, :);
        twoD_benchmark_temp = [twoD_subsample; exterior];
        score_E_benchmark_temp = input_orientedmodel(twoD_benchmark_temp, xx, yy);
        index_2 = find(score_E_benchmark_temp  > 0.9999999999999);
        twoD_benchmark = twoD_benchmark_temp(index_2, :);
        output = twoD_benchmark;
    else
        output = optimalValues;
    end
end
