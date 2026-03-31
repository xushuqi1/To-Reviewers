function results = combined_models(P, xx, yy, zz, xx1)
% COMBINED_MODELS integrates three DEA models and outputs the results of the three models
%   results = combined_models(P, xx, yy, zz, xx1)
%
% Input:
%   P   - Data matrix, column order: [X, Y, V], where X is the first group of inputs,
%         Y is desirable outputs, V is undesirable outputs.
%         X1 (the second group of inputs) is used in models 2 and 3, but the data is still
%         included in the first xx columns of P; its number of columns needs to be
%         specified separately as xx1.
%   xx  - Number of columns for the first group of inputs X (corresponds to the first xx columns of P)
%   yy  - Number of columns for desirable outputs Y (corresponds to columns xx+1 to xx+yy of P)
%   zz  - Number of columns for undesirable outputs V (corresponds to columns xx+yy+1 to xx+yy+zz of P)
%   xx1 - Number of columns for the second group of inputs X1; typically X1 is the first xx1 columns of X.
%         If omitted, defaults to xx (i.e., same as the first group of inputs).
%
% Output:
%   results - Structure containing the following fields:
%       .model13 : Efficiency values from model 1 (column vector)
%       .model15 : Efficiency values from model 2 (column vector)
%       .model16 : Efficiency values from model 3 (column vector)

    % Handle optional parameters
    if nargin < 5
        xx1 = xx;
    end

    n = size(P, 1);
    
    % ==================== model 13====================
    % Extract data
    X = P(:, 1:xx);
    Y = P(:, xx+1:xx+yy);
    V = P(:, xx+yy+1:xx+yy+zz);
    nn1 = size(X, 2);
    nn2 = size(Y, 2);
    nn3 = size(V, 2);
    
    % Preallocate
    W_m1 = zeros(3*n+1, n);
    
    for k = 1:n
        c_m1 = [zeros(3*n, 1); -1];
        A_m1 = [zeros(nn2, n), -Y', zeros(nn2, n), Y(k,:)';
                zeros(nn3, 2*n), V', V(k,:)';
                -eye(n), eye(n), zeros(n, n+1);
                -eye(n), zeros(n, n), eye(n), zeros(n, 1)];
        b_m1 = [-Y(k,:)';
                V(k,:)';
                zeros(n, 1);
                zeros(n, 1)];
        Aeq_m1 = [X', zeros(nn1, 2*n+1)];
        beq_m1 = [X(k,:)'];
        lb_m1 = zeros(3*n+1, 1);
        ub_m1 = [inf*ones(3*n, 1); 0.999999999999];
        
        % Solve using cplexlp
        W_m1(:,k) = cplexlp(c_m1, A_m1, b_m1, Aeq_m1, beq_m1, lb_m1, ub_m1);
    end
    
    L1 = W_m1(3*n+1, :);
    L_ours = 1 - L1;
    L_ours = L_ours';
    
    % ==================== model 15====================
    % Extract data (Note: X1 is the first xx1 columns of P)
    X2 = P(:, 1:xx);      % First group of inputs, used in the first linear program
    Y2 = P(:, xx+1:xx+yy);
    V2 = P(:, xx+yy+1:xx+yy+zz);
    X1_2 = P(:, 1:xx1);   % Second group of inputs, used in the second linear program
    nn1_2 = size(X2, 2);
    nn2_2 = size(Y2, 2);
    nn3_2 = size(V2, 2);
    
    % First linear program: input-oriented
    W2_m2a = zeros(n+1, n);
    for k = 1:n
        c_m2a = [zeros(n, 1); -1];
        A_m2a = [X2', zeros(nn1_2, 1);
                 -Y2', Y2(k,:)'];
        b_m2a = [X2(k,:)';
                 zeros(nn2_2, 1)];
        Aeq_m2a = [];
        beq_m2a = [];
        lb_m2a = zeros(n+1, 1);
        ub_m2a = [];
        W2_m2a(:,k) = cplexlp(c_m2a, A_m2a, b_m2a, Aeq_m2a, beq_m2a, lb_m2a, ub_m2a);
    end
    L1_m2 = W2_m2a(n+1, :);
    L1_m2 = 1 ./ L1_m2;
    
    % Second linear program: output-oriented (with undesirable output constraint)
    W2_m2b = zeros(n+1, n);
    for k = 1:n
        c_m2b = [zeros(n, 1); 1];
        A_m2b = [-X1_2', zeros(xx1, 1);
                 V2', -V2(k,:)'];
        b_m2b = [-X1_2(k,:)';
                 zeros(nn3_2, 1)];
        Aeq_m2b = [];
        beq_m2b = [];
        lb_m2b = zeros(n+1, 1);
        ub_m2b = [];
        W2_m2b(:,k) = cplexlp(c_m2b, A_m2b, b_m2b, Aeq_m2b, beq_m2b, lb_m2b, ub_m2b);
    end
    L2_m2 = W2_m2b(n+1, :);
    
    L_FGL = (L1_m2 + L2_m2) / 2;
    L_FGL = L_FGL';
    
    % ==================== model 16 ====================
    % Extract data
    X3 = P(:, 1:xx);
    X1_3 = P(:, 1:xx1);
    Y3 = P(:, xx+1:xx+yy);
    V3 = P(:, xx+yy+1:xx+yy+zz);
    
    W_m3 = zeros(2*n+1, n);
    for k = 1:n
        c_m3 = [zeros(2*n, 1); -1];
        A_m3 = [X3', zeros(xx, n+1);
                -Y3', zeros(yy, n), Y3(k,:)';
                zeros(xx1, n), -X1_3', zeros(xx1, 1);
                zeros(zz, n), V3', V3(k,:)'];
        b_m3 = [X3(k,:)';
                -Y3(k,:)';
                -X1_3(k,:)';
                V3(k,:)'];
        Aeq_m3 = [];
        beq_m3 = [];
        lb_m3 = zeros(2*n+1, 1);
        ub_m3 = [];
        W_m3(:,k) = cplexlp(c_m3, A_m3, b_m3, Aeq_m3, beq_m3, lb_m3, ub_m3);
    end
    
    L_DDF = 1 - W_m3(2*n+1, :);
    L_DDF = L_DDF';
    

    results.model13 = L_ours;
    results.model15 = L_FGL;
    results.model16 = L_DDF;
end