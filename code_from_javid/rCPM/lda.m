classdef lda < explanatory
    %LDA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mse;
        all_phenotypes;
        train_mats_weighted;
        pval;
    end
    
    methods
        
        function this = lda(group,options)
            this = this@explanatory(group,options);
            this.all_phenotypes = [];
            for i=1:length(this.phenotypes)
                this.all_phenotypes = [this.all_phenotypes this.phenotypes(i).all_behav];
            end
        end
        
        function output=run(this)
            rng(this.seed);
            x= this.all_edges;
            y = this.all_phenotypes;
            y_group = this.diagnosis;
            
            this.pval = zeros(this.num_edge,min(this.num_task,size(y,2)));
            weight = zeros(this.num_edge, this.num_task);
            
            for j_edge = 1 : this.num_edge
                edge_temp = squeeze(x(j_edge, :, 1:9));
                model= fitcdiscr(edge_temp, y_group);
                this.pval(j_edge,:) = stats.pF;
            end
            [pID,~] = this.my_fdr(this.pval(:,1),0.1);
            disp(sum(this.pval(:,1)<pID));
        end
    end
end

