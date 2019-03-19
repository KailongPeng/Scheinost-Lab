classdef group
    properties
        subjects; % N* 268*268
%         phenotypes;
        group_size;% N
        num_node;
        num_task;
        num_edge;
        all_edges;
        issym;
        k_fold; % number of folds
    end
    methods
        function this = group(subjects)
            this.subjects = subjects;
%             this.phenotypes = phenotypes;
            this.group_size = size(subjects,2);
            this.num_node = subjects(1).num_node;
            this.num_edge = subjects(1).num_edge;
            this.num_task = subjects(1).num_task;
            this.issym = subjects(1).issym;
            this.all_edges = zeros(this.num_edge,this.group_size,this.num_task);
            for i=1:this.group_size
                this.all_edges(:,i,:) = this.subjects(i).all_edges;
            end
        end
        function results = pheno_cca(this,options)
            options.subjects = this.subjects;
            for phenotype=options.phenotypes
                cca_ = cca(this,options);
                for j=1:1% FIXME change to 10
                    fprintf('Run %.0d for %s\n',j,phenotype.name) ;
                    cca_.run();
                    results(j).output = cca_.output;
                    save(sprintf('../results/RR_6tasks_%s_CCA_%.0f',phenotype.name,j),'results');
                end
            end
        end
        function results =  pheno_pca(this, options,phenotypes)
            options.subjects = this.subjects;
            options.phenotypes = this.phenotypes;
            for phenotype=phenotypes
                cpm_ = cpm(this,phenotype,options);
                for j=1:1% FIXME change to 10
                    fprintf('Run %.0d for %s\n',j,phenotype.name) ;
                    cpm_.run();
                    results(j).output = cpm_.output;
                    save(sprintf('../results/RR_6tasks_%s_PCA_%.0f',phenotype.name,j),'results');
                end
                %
            end
        end
    end
end