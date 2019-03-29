function [a, b, r, p] = mcca_eig(X, Y)
    % X should be {n_sample x n_dim1 x n_set}
    % Y should be {n_sample x n_dim2}
    n_set = size(X, 3);
    n_dx = size(X, 2);
    n_dy = size(Y, 2);
    a = zeros(n_dx, n_set);
    b = zeros(n_dy, 1);
    r = zeros(n_set, 1);
    p = zeros(n_set, 1);
    
    % precalculate necessary variables
    X_center = X - mean(X, 1);
    Y_center = Y - mean(Y, 1);
    rab = zeros(n_dx, n_dy, n_set);
    raa = zeros(n_dx, n_dx, n_set);
    rbb = (Y-mean(Y))' * (Y-mean(Y));
    rbb_inv = pinv(rbb);
    for i_set = 1 : n_set
        temp_X = X_center(:, :, i_set);
        rab(:, :, i_set) = temp_X'*Y_center;
        raa(:, :, i_set) = temp_X'*temp_X;
    end 
    
    S = spalloc(n_dx*n_set+n_dy, n_dx*n_set+n_dy, 2*n_set*n_dx*n_dy);
    for i_set = 1 : n_set
        S(n_set*n_dx+1:end, (i_set-1)*n_dx+1:i_set*n_dx) = rbb_inv*rab(:, :, i_set)';
        S((i_set-1)*n_dx+1:i_set*n_dx, n_set*n_dx+1:end) = pinv(raa(:, :, i_set))*rab(:, :, i_set);
    end
    [V, ~] = eigs(S, 1);

%     S = spalloc(n_dx*n_set+n_dy, n_dx*n_set+n_dy, 2*n_set*n_dx*n_dy);
%     D = spalloc(n_dx*n_set+n_dy, n_dx*n_set+n_dy, n_dx*n_dx*n_set+n_dy*n_dy);
%     for i_set = 1 : n_set
%         S(n_set*n_dx+1:end, (i_set-1)*n_dx+1:i_set*n_dx) = rab(:, :, i_set)';
%         S((i_set-1)*n_dx+1:i_set*n_dx, n_set*n_dx+1:end) = rab(:, :, i_set);
%         D((i_set-1)*n_dx+1:i_set*n_dx, (i_set-1)*n_dx+1:i_set*n_dx) = raa(:, :, i_set);
%     end
%     D(n_set*n_dx+1:end, n_set*n_dx+1:end) = rbb;
%     [V, ~] = eigs(S, D, 1);
    
    a = reshape(V(1:n_dx*n_set), n_dx, n_set);
    b = V(n_dx*n_set+1:end);
    for i_set = 1 : n_set
%         temp_numr = (X_center(:, :, i_set)*a(:, i_set))'*(Y_center*b);
%         temp_denom = sqrt(a(:, i_set)'*X_center(:, :, i_set)'*X_center(:, :, i_set)*a(:, i_set))*sqrt((Y_center*b)'*Y_center*b);
%         r(i_set) = temp_numr/temp_denom;
        [r(i_set), p(i_set)] = corr((X_center(:, :, i_set)*a(:, i_set)), (Y_center*b));
    end
end