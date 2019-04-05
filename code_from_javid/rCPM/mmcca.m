classdef mmcca < predictory
    %mmcca Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num_behav;
        num_factors;
        thresh1;
        thresh2;
    end
    
    methods
        function this = mmcca(group,phenotype,options)
            %mmcca Construct an instance of this class
            %   Detailed explanation goes here
            this = this@predictory(group,phenotype,options);
            this.all_phenotypes = [];
            for i=1:length(this.phenotypes)
                this.all_phenotypes = [this.all_phenotypes this.phenotypes(i).all_behav];
            end
            this.num_behav = size(this.phenotype.all_behav,2);
            
            if isfield(options,'threshold1')
                this.threshold1=options.threshold1;
            else
                this.threshold1=0.3;
            end
            if isfield(options,'threshold2')
                this.threshold2=options.threshold2;
            else
                this.threshold2=0.1;
            end
        end
        function output=run(this)
            y = zeros(this.num_sub_total, 1);
            this.new_behav = zeros(this.num_sub_total, 1);
            rng(this.seed);
            indices = crossvalind('Kfold', this.num_sub_total, k);
            this.r_pearson = zeros(k, 1);
            this.r_rank = zeros(k, 1);
            this.mse = zeros(k, 1);
            this.q_s = zeros(k, 1);
            this.all_edge_weight = zeros(k, this.num_edge);
            this.all_task_weight = zeros(k, this.num_task);
            this.all_behav_weight = zeros(k, this.num_behav);
            this.lambda_total = zeros(1, k); % store all the lambda 
            %%%%%%%%%%%% make sure we have not initialized in predictory's
            %%%%%%%%%%%% constructor
            
            for i_fold = 1 : k
                tStart = tic;
                fprintf('%dth fold\n', i_fold);
 
                test.indx = (indices == i_fold);
                train.indx = ~(indices == i_fold);
                test.x = this.all_edges(:,test.indx);
                train.x = this.all_edges(:,train.indx);
                test.y = this.phenotype.all_behav(indices == i_fold,:);
                train.y = this.phenotype.all_behav(~(indices == i_fold),:);
                
                num_train = size(train.y, 1); % FIXME double check
                num_test = size(test.x, 2); % FIXME double check
                
                % select edges that are significant in every task
                all_p = zeros(this.num_edge, this.num_behav, this.num_task);
                for j_task = 1 : num_task
                    [~, all_p(:, :, j_task)] = corr(train_mats(:, :, j_task)', train_behav);
                end
                all_p = mean(mean(all_p, 2),3);
                edge_idx1 = find(all_p<this.thresh1);
                train_mats = train_mats(edge_idx1, :, :);
                
                disp(numel(edge_idx1))
                
                [task_weight, behav_weight, ~, edge_p] = mcca_eig(permute(train_mats, [2, 3, 1]), train_behav);                
                
                mean_feature = squeeze(mean(train_mats, 2));
                train_behav = (train_behav-mean(train_behav, 1)) * behav_weight;
                test_behav = (test_behav-mean(test_behav, 1)) * behav_weight;
                
                % get sum of all edges in TRAIN subs
                all_edges_new = zeros(numel(edge_idx1), this.num_sub_total);
                for j_sub = 1 : this.num_sub_total
                    all_edges_new(:, j_sub) = sum((squeeze(this.all_edges(edge_idx1, j_sub, :)) - mean_feature) .* task_weight', 2);
                end
                train_mats = all_edges_new(:, train_idx);
                test_mats = all_edges_new(:, test_idx);
                
                % define masks
                edge_idx2 = edge_p < this.thresh2;
                
                if ~exist('this.lambda', 'var')  || length(this.lambda) ~= 1
                    [fit_coef, fit_info] = lasso(train_mats(edge_idx2, :)', train_behav, 'Alpha',1e-6, 'CV', 10);
                    idxLambda1SE = fit_info.Index1SE;
                    coef = fit_coef(:,idxLambda1SE);
                    coef0 = fit_info.Intercept(idxLambda1SE);
                    this.lambda_total(i_fold) = fit_info.Lambda(idxLambda1SE);
                else
                    [coef, fit_info] = lasso(train_mats(edge_idx2, :)', train_behav, 'Alpha',1e-6, 'Lambda', lambda);
                    coef0 = fit_info.Intercept;
                end
                

            end
        end
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

            a = reshape(V(1:n_dx*n_set), n_dx, n_set);
            b = V(n_dx*n_set+1:end);
            for i_set = 1 : n_set
                [r(i_set), p(i_set)] = corr((X_center(:, :, i_set)*a(:, i_set)), (Y_center*b));
            end
        end
        function evaluate()
                            this.y(test_idx) = test_mats(edge_idx2, :)'*coef+coef0;
                new_behav(test_idx) = test_behav;
                % compare predicted and observed behaviors
                [r_pearson(i_fold), ~] = corr(this.y(test_idx), test_behav);
                [r_rank(i_fold), ~] = corr(this.y(test_idx), test_behav, 'type', 'spearman');
                this.mse(i_fold) = sum((this.y(test_idx) - test_behav).^2) / num_test;
                this.q_s(i_fold) = 1 - this.mse(i_fold) / var(test_behav, 1);
                
                
                % calculate or save out some weights
                all_std = std(permute(all_edges(:, train_idx, :), [1,3,2]),[], 3);
                all_std = all_std(edge_idx1(edge_idx2), :);
                this.all_task_weight(i_fold, :) =  abs(diag(task_weight(:, edge_idx2)*(all_std.*abs(coef))));
                
                this.all_edge_weight(i_fold, edge_idx1(edge_idx2)) = std(train_mats(edge_idx2, :), [], 2).*abs(coef);
                
                this.all_behav_weight(i_fold, :) = behav_weight;
                
                tElapsed = toc(tStart)
        end
    end
end

