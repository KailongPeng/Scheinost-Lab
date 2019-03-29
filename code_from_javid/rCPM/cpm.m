classdef cpm < cbase
    methods
        function this = cpm(group,phenotype,options)
            this = this@cbase(group,phenotype,options);
        end
        function run(this)
            
            all_edges = this.group.all_edges;
            all_edges = permute(all_edges, [1, 3, 2]);
            all_edges = reshape(all_edges, [],this.num_sub_total);
            
            y = zeros(this.num_sub_total, 1);
            rng(this.seed);
            indices = crossvalind('Kfold', this.num_sub_total,this.k);
            
            for i_fold = 1 : this.k
                fprintf('%dth fold\n', i_fold);
                test.indx = (indices == i_fold);
                train.indx = ~(indices == i_fold);
                test.x = this.all_edges(:,test.indx);
                train.x = this.all_edges(:,train.indx);
                test.y = this.phenotype.all_behav(indices == i_fold,:);
                train.y = this.phenotype.all_behav(~(indices == i_fold),:);
                
                % first step univariate edge selection
                [~, edge_p] = corr(train.x', train.y);
                edges_1 = find(edge_p < this.thresh);
                
                % build model on TRAIN subs
                if ~exist('lambda', 'var')  || length(this.lambda) ~= 1
                    [fit_coef, fit_info] = lasso(train.x(edges_1, :)', train.y, 'Alpha',this.alpha, 'CV', 10);
                    idxLambda1SE = fit_info.Index1SE;
                    coef = fit_coef(:,idxLambda1SE);
                    coef0 = fit_info.Intercept(idxLambda1SE);
                    this.lambda_total(i_fold) = fit_info.Lambda(idxLambda1SE);
                else
                    [coef, fit_info] = lasso(train.x(edges_1, :)', train.y, 'Alpha',this.alpha, 'Lambda', this.lambda);
                    coef0 = fit_info.Intercept;
                end
                
                % run model on TEST sub with the best lambda parameter
                test.x = all_edges(:, test.indx);
                y(test.indx) = test.x(edges_1, :)'*coef+coef0;
                
                this.coef_total(edges_1, i_fold) = coef;
                this.coef0_total(:, i_fold) = coef0;
                
            end
            
            % compare predicted and observed behaviors
            [r_pearson, ~] = corr(y, this.phenotype.all_behav);
            [r_rank, ~] = corr(y, this.phenotype.all_behav, 'type', 'spearman');
            mse = sum((y - this.phenotype.all_behav).^2) / this.group.group_size;
            q_s = 1 - mse / var(this.phenotype.all_behav, 1);
        end
        
        function run_diagnosis(this)
            rng(this.seed);
            indices = crossvalind('Kfold', this.num_sub_total, this.k);
            for i_fold = 1 : this.k
                fprintf('%dth fold\n', i_fold);
                test.indx = (indices == i_fold);
                train.indx = ~(indices == i_fold);
                test.x = this.all_edges(:,test.indx);
                train.x = this.all_edges(:,train.indx);
                test.y = this.phenotype.all_behav(indices == i_fold,:);
                train.y = this.phenotype.all_behav(~(indices == i_fold),:);
                
                [train_behav,test_behav] = this.kfolds_pca(train,test,i_fold);
                this.phenotype.pca_behav(test.indx,:) = test_behav;
                % first step univariate edge selection
                [~, edge_p] = corr(train.x', train_behav);
                edges_1 = find(edge_p < this.thresh);
                % build model on TRAIN subs
                if ~exist('lambda', 'var')  || length(this.lambda) ~= 1
                    [fit_coef, fit_info] = lasso(train.x(edges_1, :)', train_behav, 'Alpha',this.v_alpha, 'CV', this.k);
                    idxLambda1SE = fit_info.Index1SE;
                    coef = fit_coef(:,idxLambda1SE);
                    coef0 = fit_info.Intercept(idxLambda1SE);
                    this.lambda_total(i_fold) = fit_info.Lambda(idxLambda1SE);
                else
                    [coef, fit_info] = lasso(train.x(edges_1, :)', train_behav, 'Alpha',this.v_alpha, 'Lambda', this.lambda);
                    coef0 = fit_info.Intercept;
                end
                % run model on TEST sub with the best lambda parameter
                this.Y(test.indx) = test.x(edges_1, :)'*coef+coef0; % this might need to be changed
                this.coef_total(edges_1, i_fold) = coef;
                this.coef0_total(:, i_fold) = coef0;
            end
        end
        function evaluate(this)
            [this.r_pearson, ~] = corr(this.Y, this.phenotype.pca_behav);
            [this.r_rank, ~] = corr(this.Y, this.phenotype.pca_behav, 'type', 'spearman');
            this.mse = sum((this.Y - this.phenotype.pca_behav).^2) / this.num_sub_total;
            this.q_s = 1 - this.mse / var(this.Y, 1);
            fprintf('q_s=%f\n',this.q_s);
        end
    end
end