classdef subject < handle
    properties
        
        bart; % 268*268 connectome of resting
        pamenc; % 268*268 connectome of resting
        pamret; % 268*268 connectome of working memory
        scap; % 268*268 connectome of gambling
        stopsignal; % 268*268 connectome of motor
        taskswitch; % 268*268 connectome of language
        
        rest1; % 268*268 connectome of resting
        rest2; % 268*268 connectome of resting
        wm; % 268*268 connectome of working memory
        gambling; % 268*268 connectome of gambling
        motor; % 268*268 connectome of motor
        language; % 268*268 connectome of language
        social; % 268*268 connectome of social cognition
        relational; % 268*268 connectome of relational processing
        emotion; % 268*268 connectome of emotion
        
        all_edges; % 6 * 268*268
        Matrix; % 268*268 matrix based on current task
        gender; % male: 1, female: 0
        age; % 
        dim; % 268*268
        
        id;
        issym;
        num_node;
        num_task;
        num_edge;
    end
    methods
        function this = subject(x,id,dataset, gender,age,mask)
            if nargin > 1
                this.id = id;
                this.num_node = size(x, 1);
                this.num_task = size(x, 4);
                is_sym = issymmetric(x(:, :, 1, 1));
                if is_sym
                    this.num_edge = this.num_node * (this.num_node - 1) / 2;
                else
                    this.num_edge = this.num_node * this.num_node;
                end
                % apply mask if it has been set
                x = squeeze(x);
                if(dataset=="UCLA")
                    this.dim = size(x(:,:,1));
                    this.bart = reshape(x(:,:,1),[],1);
                    this.pamenc = reshape(x(:,:,2),[],1);
                    this.pamret = reshape(x(:,:,3),[],1);
                    this.scap = reshape(x(:,:,4),[],1);
                    this.stopsignal = reshape(x(:,:,5),[],1);
                    this.taskswitch = reshape(x(:,:,6),[],1);
                elseif(dataset=="HCP")
                    this.dim = size(x(:,:,1));
                    this.gambling = reshape(x(:,:,1),[],1);
                    this.rest1 = reshape(x(:,:,2),[],1);
                    this.rest2 = reshape(x(:,:,3),[],1);
                    this.language = reshape(x(:,:,4),[],1);
                    this.motor = reshape(x(:,:,5),[],1);
                    this.relational = reshape(x(:,:,6),[],1);
                    this.social = reshape(x(:,:,7),[],1);
                    this.wm = reshape(x(:,:,8),[],1);
                    this.emotion = reshape(x(:,:,9),[],1);
                    this.gender = gender;
                    this.age = age;
                end
                this.all_edges = zeros(this.num_edge,this.num_task);
                    for j_task = 1 : this.num_task
                        if is_sym
                            this.all_edges(:, j_task) = squareform(tril(x(:, :, j_task), -1));
                        else
                            this.all_edges(:, j_task) = reshape(x(:, :, j_task), [], 1);
                        end
                    end
%                 this.all_edges = permute(this.all_edges, [1, 3, 2]);
%                 this.all_edges = reshape(this.all_edges, [], 1);
                this.issym = issymmetric(x(:, :, 1, 1));
                this.applyMask(mask);
            end
        end
        function this=applyMask(this,mask)
            if(mask)
                missing_nodes = [60, 100, 108, 109, 112, 115,116, 118, 129, 189, ...
                    202, 239, 240, 242, 243, 249, 250, 266];
                f1 = load('../data/abby_all_task_edges_pos_p0.001.mat');
                f2 = load('../data/abby_all_task_edges_neg_p0.001.mat');
                mask = squeeze(f1.task_pos_edge(2,:,:)+f2.task_neg_edge(2,:,:)); % wm mask
                for i=missing_nodes
                    b1 = zeros(1,size(mask,1));
                    b2 = zeros(size(mask,1)+1,1);
                    mask = [mask(1:i,:); b1; mask(i+1:end,:)];
                    mask = [mask(:,1:i) b2 mask(:,i+1:end)];                 
                end % now we have 268*268

                this.bart = reshape(mask,[],1).*this.bart;
                this.pamenc = reshape(mask,[],1).*this.pamenc;
                this.pamret = reshape(mask,[],1).*this.pamret;
                this.scap = reshape(mask,[],1).*this.scap;
                this.stopsignal = reshape(mask,[],1).*this.stopsignal;
                this.taskswitch = reshape(mask,[],1).*this.taskswitch;
            end
        end
    end
end