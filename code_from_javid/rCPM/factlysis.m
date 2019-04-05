classdef factlysis < mmcca
    %FACTLYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    % properties
    %     Property1
    % end
    
    methods
        function this = factlysis(group,phenotype,options)
            %FACTLYSIS Construct an instance of this class
            %   Detailed explanation goes here
            this = this@mmcca(group,phenotype,options);
            
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
    end
end

