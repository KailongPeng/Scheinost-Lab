classdef cca < predictory
    properties
        mse;
        all_phenotypes;
        train_mats_weighted;
        pval;
    end
    methods
        function this = cca(group,options)
            this = this@predictory(group,options);
            this.all_phenotypes = [];
            for i=1:length(this.phenotypes)
                this.all_phenotypes = [this.all_phenotypes this.phenotypes(i).all_behav];
            end
        end
        
        function output=run(this)
            rng(this.seed);
%             indices = crossvalind('Kfold', this.num_sub_total, this.k);
%             mask = zeros(this.num_node, this.num_node, this.k); %store all the edges selected
            x= this.all_edges;
            y = this.all_phenotypes;
            y_group = this.diagnosis; 
            
            % weight calculation
            this.pval = zeros(this.num_edge,min(this.num_task,size(y,2)));
            weight = zeros(this.num_edge, this.num_task);
            
            for j_edge = 1 : this.num_edge
                %disp(j/num_edge)
                edge_temp = squeeze(x(j_edge, :, :));
%                 if rank(edge_temp)==this.num_task
                    [A, B, ~, ~, ~, stats] = canoncorr(edge_temp, y);
%                     if B > 0
%                         weight(j_edge, :) = A;
%                     elseif B<0
%                         weight(j_edge, :) = -A;
%                     end
                    this.pval(j_edge,:) = stats.pF;
%                 end
            end
            % define masks
            [pID,~] = this.my_fdr(this.pval(:,1),0.1);
            disp(sum(this.pval(:,1)<pID));
%             edge_idx = edge_p < P;
%             output.edge_ix = edge_idx;
%             edge_idx = edge_idx(:,1); % only on first dimension
        end
        
    end
end