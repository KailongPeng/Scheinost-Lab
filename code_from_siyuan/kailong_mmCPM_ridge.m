function [q_s, r_pearson, r_rank, y, new_behav, all_edge_weight, all_behav_weight, all_task_weight, lambda_total] = mmCPM_ridge(all_mats, all_behav, thresh1, thresh2, lambda, k, seed)
    %CPM Connectome-based predictive modeling using univariate feature selection 
    %
    %   [q_s, r_pearson, r_rank, y, mask] = mCPM(all_mats, all_behav, 0.1)
    %   [q_s, r_pearson, r_rank, y, mask] = mCPM(all_mats, all_behav, 0.1, 10, 665)
    %
    %   Input:      all_mats,           connectome of all the subjects and tasks
    %                                   [regions x regions x subjects x tasks]
    %
    %               all_behav,          behavior of all the subjects
    %                                   [subjects x behavior]
    %
    %               thresh,             p-value threshold used for univariate
    %                                   feature selection
    %
    %               k(optional),        number of folds in k-fold cross
    %                                   validation, default is 10
    %
    %               seed(optional),     random seed, default is 665
    %
    %   Output:     q_s,                cross-validated R^2 between predicted
    %                                   value with ground truth (all_behav)
    %
    %               r_pearson,          cross-validated pearson correlation 
    %                                   between predicted value with ground 
    %                                   truth (all_behav)
    %
    %               r_rank,             cross-validated spearman correlation
    %                                   between predicted value with ground 
    %                                   truth (all_behav)
    %
    %               y,                  predicted value for all the subjects
    % 
    %               mask,               network selected
    %
    %   Siyuan Gao, Yale University, 2019-2020
    
    %% initialization
    if nargin < 3
        error('not enough arguments, please check the help')
    end
    
    if ~exist('k', 'var')
        k = 10;
    end
    
    if ~exist('seed', 'var')
        seed = 665;
    end    
    
    num_sub_total = size(all_mats, 3);
    num_node = size(all_mats, 1);
    num_task = size(all_mats, 4);
    num_behav = size(all_behav, 2);
    is_sym = issymmetric(all_mats(:, :, 1, 1));
    if is_sym
        num_edge = num_node * (num_node - 1) / 2;
    else
        num_edge = num_node * num_node;
    end
    
    %% convert connectivity to edge matrix
    all_edges = zeros(num_edge, num_sub_total, num_task);
    for i_sub = 1 : num_sub_total
        for j_task = 1 : num_task
            if is_sym
                all_edges(:, i_sub, j_task) = squareform(tril(all_mats(:, :, i_sub, j_task), -1));
            else
                all_edges(:, i_sub, j_task) = reshape(all_mats(:, :, i_sub, j_task), [], 1);
            end
        end
    end

    
    %% main
    num_edge = size(all_edges, 1);
    y = zeros(num_sub_total, 1);
    new_behav = zeros(num_sub_total, 1);
    rng(seed);
    indices = crossvalind('Kfold', num_sub_total, k);
    r_pearson = zeros(k, 1);
    r_rank = zeros(k, 1);
    mse = zeros(k, 1);
    q_s = zeros(k, 1);
    all_edge_weight = zeros(k, num_edge);
    all_task_weight = zeros(k, num_task);
    all_behav_weight = zeros(k, num_behav);
    lambda_total = zeros(1, k); % store all the lambda
    for i_fold = 1 : k
        tStart = tic;
        fprintf('%dth fold\n', i_fold);
        
        test_idx = (indices==i_fold);
        train_idx = (indices~=i_fold);
        train_mats = all_edges(:, train_idx, :);
        train_behav = all_behav(train_idx, :);
        test_behav = all_behav(test_idx, :);
        test_mats = all_edges(:, test_idx, :);
        num_train = size(train_behav, 1);
        num_test = size(test_mats, 2);
        
        % pca on train data, get coefficient; when test model, apply
        % coeffient to test data to test
        [coeff,score,latent,tsquared,explained,mu] = pca(train_behav,'algorithm','als');
        train_behav = score;
        % factor analysis on train data, get coefficient; when test model, apply
        % coeffient to test data to test
        
                
        % select edges that are significant in every task
        all_p = zeros(num_edge, num_behav, num_task);
        for j_task = 1 : num_task
            [~, all_p(:, :, j_task)] = corr(train_mats(:, :, j_task)', train_behav);
        end
        %         all_p = max(max(all_p, [], 2), [], 3);
        all_p = mean(mean(all_p, 2),3);
        edge_idx1 = find(all_p<thresh1);
        train_mats = train_mats(edge_idx1, :, :);
        
        disp(numel(edge_idx1))
        
        [task_weight, behav_weight, ~, edge_p] = mcca_eig(permute(train_mats, [2, 3, 1]), train_behav);
%         [task_weight, behav_weight, r, edge_p, a_list, b_list, r_list, stop_iter] = mcca_iter(permute(train_mats, [2, 3, 1]), train_behav, 1000, 1e-5);
        
        
        mean_feature = squeeze(mean(train_mats, 2));
        train_behav = (train_behav-mean(train_behav, 1)) * behav_weight;
        test_behav = (test_behav-mean(test_behav, 1)) * behav_weight;     

        % get sum of all edges in TRAIN subs
        all_edges_new = zeros(numel(edge_idx1), num_sub_total);
        for j_sub = 1 : num_sub_total
            all_edges_new(:, j_sub) = sum((squeeze(all_edges(edge_idx1, j_sub, :)) - mean_feature) .* task_weight', 2);
        end
        train_mats = all_edges_new(:, train_idx);
        test_mats = all_edges_new(:, test_idx);
        
        % define masks
        edge_idx2 = edge_p < thresh2;
        
        if ~exist('lambda', 'var')  || length(lambda) ~= 1
            [fit_coef, fit_info] = lasso(train_mats(edge_idx2, :)', train_behav, 'Alpha',1e-6, 'CV', 10);
            idxLambda1SE = fit_info.Index1SE;
            coef = fit_coef(:,idxLambda1SE);
            coef0 = fit_info.Intercept(idxLambda1SE);
            lambda_total(i_fold) = fit_info.Lambda(idxLambda1SE);
        else
            [coef, fit_info] = lasso(train_mats(edge_idx2, :)', train_behav, 'Alpha',1e-6, 'Lambda', lambda);
            coef0 = fit_info.Intercept;
        end
        y(test_idx) = test_mats(edge_idx2, :)'*coef+coef0;
        new_behav(test_idx) = test_behav;
        % compare predicted and observed behaviors
        [r_pearson(i_fold), ~] = corr(y(test_idx), test_behav);
        [r_rank(i_fold), ~] = corr(y(test_idx), test_behav, 'type', 'spearman');
        mse(i_fold) = sum((y(test_idx) - test_behav).^2) / num_test;
        q_s(i_fold) = 1 - mse(i_fold) / var(test_behav, 1);
        
        
        % calculate or save out some weights
        all_std = std(permute(all_edges(:, train_idx, :), [1,3,2]),[], 3);
        all_std = all_std(edge_idx1(edge_idx2), :);
        all_task_weight(i_fold, :) =  abs(diag(task_weight(:, edge_idx2)*(all_std.*abs(coef))));
        
        all_edge_weight(i_fold, edge_idx1(edge_idx2)) = std(train_mats(edge_idx2, :), [], 2).*abs(coef);
        
        all_behav_weight(i_fold, :) = behav_weight;
        
        tElapsed = toc(tStart)
    end
    
end

