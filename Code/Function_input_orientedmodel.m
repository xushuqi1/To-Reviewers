function [output] = Function_model(P, xx, yy)
    X = P(:, 1:xx);
    Y = P(:, xx+1:xx+yy);  
    n = size(X, 1);

    optimalValues = zeros(1, n);

    parfor k = 1:n
        c = [zeros(n, 1); 1];
        A = [X', -X(k, :)'; -Y', zeros(yy, 1)];
        b = [zeros(xx, 1); -Y(k, :)'];
        Aeq = [ones(1, n), 0];
        beq = 1;
        lb = zeros(n+1, 1);
        ub = [];
        op = optimoptions('linprog', 'display', 'none');
        
        [sol, fval] = linprog(c, A, b, Aeq, beq, lb, ub, op);

        optimalValues(k) = fval;
    end

    output = optimalValues;
end
