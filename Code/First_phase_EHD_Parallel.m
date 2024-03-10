function [output] = First_phase_EHD_Parallel(P, xx, yy, a)        %  P refers to a n*(xx+yy) matrix.
% a=int(n^0.5)
% The input-oriented model is used here.
% To avoid the effect of the precision setting, the linear programming we use in determining the interior
% points comes from Dula (2011)
    E_EHD = EHD2019_initial_step(P, xx, yy, a);            
    score_E_EHD = input_orientedmodel(E_EHD, xx, yy);
    index_1 = find(score_E_EHD > 0.9999999999999);
    EHD_subsample = E_EHD(index_1, :);

    n = size(EHD_subsample, 1);    
    nnn = size(P, 1);
    X1 = EHD_subsample(:, 1:xx);
    Y1 = EHD_subsample(:, xx+1:xx+yy);

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
        EHD_benchmark_temp = [EHD_subsample; exterior];
        score_E_benchmark_temp = input_orientedmodel(EHD_benchmark_temp, xx, yy);
        index_2 = find(score_E_benchmark_temp > 0.9999999999999);
        EHD_benchmark = EHD_benchmark_temp(index_2, :);
        output = EHD_benchmark;
    else
        output = optimalValues;
    end
end
