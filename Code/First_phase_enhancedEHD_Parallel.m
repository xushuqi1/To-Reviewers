function [output] = First_phase_enhancedEHD_Parallel(P, xx, yy, a)              % a=int(n^0.5) 

    E_EHD = Function_initial_step_enhanced_EHD (P, xx, yy, a);            
    score_E_EHD = input_orientedmodel(E_EHD, xx, yy);
    index_1 = find(score_E_EHD > 0.9999999999999);
    EHD_subsample = E_EHD(index_1, :);

    n = size(EHD_subsample, 1);    
    nnn = size(P, 1);
    X1 = EHD_subsample(:, 1:xx);
    Y1 = EHD_subsample(:, xx+1:xx+yy);

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
        EHD_benchmark_temp = [EHD_subsample; exterior];
        score_E_benchmark_temp = input_orientedmodel(EHD_benchmark_temp, xx, yy);
        index_2 = find(score_E_benchmark_temp > 0.9999999999999);
        EHD_benchmark = EHD_benchmark_temp(index_2, :);
        output = EHD_benchmark;
    else
        output = optimalValues;
    end
end
