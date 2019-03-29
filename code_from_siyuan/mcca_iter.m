function [a, b, r, p, a_list, b_list, r_list, stop_iter] = mcca_iter(X, Y, n_iter, eps)
    if ~exist('eps', 'var')
        eps = 1e-5;
    end
    % X should be {n_sample x n_dim1 x n_set}
    % Y should be {n_sample x n_dim2}
    n_set = size(X, 3);
    n_dx = size(X, 2);
    n_dy = size(Y, 2);
    a = zeros(n_dx, n_set);
    b = zeros(n_dy, 1);
    r = zeros(n_set, 1);
    p = zeros(n_set, 1);
    for i_set = 1 : n_set
        [aa, bb] = canoncorr(X(:, :, i_set), Y);
        a(:, i_set) = aa(:, 1);
        b = bb(:, 1);
    end
    % precalculate necessary variables
    X_center = X - mean(X, 1);
    Y_center = Y - mean(Y, 1);
    rab = zeros(n_dx, n_dy, n_set);
    raa = zeros(n_dx, n_dx, n_set);
    rbb = (Y-mean(Y))' * (Y-mean(Y));
    for i_set = 1 : n_set
        temp_X = X_center(:, :, i_set);
        rab(:, :, i_set) = temp_X'*Y_center;
        raa(:, :, i_set) = temp_X'*temp_X;
    end    
    
    
    % interleaving update the parameters, to find the fix point
    a_list = zeros(n_iter+1, n_dx, n_set);
    b_list = zeros(n_iter+1, n_dy);
    r_list = zeros(n_iter+1, n_set);
    
    a_list(1, :, :) = a;
    b_list(1, :) = b;
    for i_set = 1 : n_set
        r_list(1, i_set) = corr((X_center(:, :, i_set)*a(:, i_set)), (Y_center*b));
    end
    
    for iter = 1 : n_iter
%         disp(iter)
        %calculate necessary variable
        temp_sum = 0;
        for i_set = 1 : n_set
             temp_sum = temp_sum + rab(:, :, i_set)'*a(:, i_set);
        end
        %update parameter
        b = pinv(rbb)*temp_sum;
        for i_set = 1 : n_set
            a(:, i_set) = pinv(raa(:, :, i_set))*rab(:, :, i_set)*b;
        end
        scale = 0;
        for i_set = 1 : n_set
            scale = a(:, i_set)'*raa(:, :, i_set)*a(:, i_set);
        end
        scale = scale + b'*rbb*b;
        a = a./sqrt(scale);
        b = b./sqrt(scale);

        a_list(iter+1, :, :) = a;
        b_list(iter+1, :) = b;
        for i_set = 1 : n_set
            r_list(iter+1, i_set) = corr((X_center(:, :, i_set)*a(:, i_set)), (Y_center*b));
        end
%         if iter>1 && norm(a-squeeze(aaa(iter-1, :, :)),'fro') < eps
        if iter>1 && abs(1-corr(b, b_list(iter-1, :)')) < eps
            disp('stopped early')
            break
        end
        
    end
    for i_set = 1 : n_set
        [r(i_set), p(i_set)] = corr((X_center(:, :, i_set)*a(:, i_set)), (Y_center*b));
    end
    stop_iter = iter;
end