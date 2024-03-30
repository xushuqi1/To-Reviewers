function [output] = Function_First_phase_twoD_Parallel(P, xx, yy) %  P refers to a n*(xx+yy) matrix.
% The input-oriented model is used here.  

    
    twoD_subsample = Function_twoD(P, xx, yy);
    n = size(twoD_subsample, 1);      
    nnn = size(P, 1);
    X1 = twoD_subsample(:, 1:xx);
    Y1 = twoD_subsample(:, xx+1:xx+yy);

    unsolvableOrPositiveK = [];
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
                unsolvableOrPositiveK = [unsolvableOrPositiveK, k];
            end
            optimalValues(k) = fval;
        catch
            unsolvableOrPositiveK = [unsolvableOrPositiveK, k];
        end
    end

    if ~isempty(unsolvableOrPositiveK)
        unsolvableOrPositiveK = unsolvableOrPositiveK';
        exterior = P(unsolvableOrPositiveK, :);
        twoD_benchmark_temp = [twoD_subsample; exterior];
        score_E_benchmark_temp = input_orientedmodel(twoD_benchmark_temp, xx, yy);
        index_2 = find(score_E_benchmark_temp  > 0.9999999999999);
        twoD_benchmark = twoD_benchmark_temp(index_2, :);
        output = twoD_benchmark;
    else
        output = optimalValues;
    end
end
