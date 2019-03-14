classdef manova < explanatory
    %MANOVA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dimensions;
        pval;
        statistics;
        
        
    end
    
    methods
        function this = manova(group,options)
            %MANOVA Construct an instance of this class
            %   Detailed explanation goes here
            addpath ../lib;
            this = this@explanatory(group,options);
        end
        function output=HottelingT(this)
            addpath ../HotellingT2;
            
            inp2 = this.diagnosis;
            male = find(inp2>0);
            female = find(inp2==0);
            n1= size(male,1);
            n2=size(female,1);
            this.pval = zeros(this.num_edge,length(unique(inp2))-1);
            posu = ones(this.num_edge,1);
            negu = ones(this.num_edge,1);
            pval_all = zeros(this.num_edge,1);
            

            for j_edge = 1 : this.num_edge
                
                a=squeeze(this.all_edges(j_edge, :,:));
                %                 a(:,9)=[];
                inp1=a;
                %                 inp1 = squeeze(this.all_edges(j_edge, :,:));
                X = [inp1(male,:);inp1(female,:)];
                [p] = T2Hot2ihe(X,n1,n2,0.05);
                [h,pvalu,ci,stats]=ttest2(inp1(male,:),inp1(female,:));
             %   if(pvalu<0.05 && stats.tstat>0)
              %      posu(j_edge)=pvalu;
               % elseif(pvalu<0.05 && stats.tstat<0)
                %    negu(j_edge)=pvalu;
                %end
                this.pval(j_edge,:) = p(1);
                %pval_all(j_edge,:)=pvalu;
            end
            [pID,~] = this.my_fdr(this.pval(:,1),0.005);
            [pID_ttest,~] = this.my_fdr(pval_all,0.005);
            
            [N_CNT,CON_MAT,PVAL]=this.myNBSfdr(this.pval(:,1),100);
            NBS_CORRECTED= full(cell2mat(CON_MAT));
            
            network = zeros(268,268);
            pos = zeros(268,268);
            neg = zeros(268,268);
            
            lindices = tril(ones(268));
            lindices(find(eye(268)))=0;
            lindices = find(lindices);
            
            uindices = triu(ones(268));
            uindices(find(eye(268)))=0;
            uindices = find(uindices);
            significants = this.pval(:,1)<pID;
            
            network(lindices) = significants;
            network = network + network';
%             network(uindices) = significants;
            
            pos(lindices) = posu<pID_ttest;
            pos=pos+pos';
%             pos(uindices) = posu<pID_ttest;
            neg(lindices) = negu<pID_ttest;
            neg = neg+neg';
%             neg(uindices) = negu<pID_ttest;
            
%             csvwrite('/Users/javid/Dropbox/PhD/visualization/server/data/task+rest1/sample.gam.csv', [1:268;network]);
            csvwrite('/Users/javid/Deesktop/pos.1.csv', pos);
            csvwrite('/Users/javid/Desktop/neg.1.csv', neg);
            
            disp(sum(this.pval(:,1)<pID));
            disp(sum(pos,'all')/2);
            disp(sum(neg,'all')/2);
            output = this.dimensions;
        end
        function output=run(this)
            inp2 = this.diagnosis;
            this.pval = zeros(this.num_edge,length(unique(inp2))-1);
            for j_edge = 1 : this.num_edge
                inp1 = squeeze(this.all_edges(j_edge, :,6))';
                [d,p,~] = manova1(inp1,inp2);
                this.dimensions(j_edge) = d;
                this.pval(j_edge,:) = p(1);
            end
            [pID,pN] = this.my_fdr(this.pval(:,1),0.01);
            [N_CNT,CON_MAT,PVAL]=this.myNBSfdr(this.pval(:,1),20);
            NBS_CORRECTED= full(cell2mat(CON_MAT));
            output = this.dimensions;
        end
        function [N_CNT,CON_MAT,PVAL]=myNBSfdr(this,network,iterations)
            [matrix,biggestArea,lcc] = this.permTest(this.pval(:,1),iterations);
            stats.alpha = 0.005;
             stats.N=this.num_node;
             stats.test_stat=matrix;
             [N_CNT,CON_MAT,PVAL]=NBSfdr(stats);
%             stats.size=;
            
        end

        function [matrix,biggestArea,lcc]=permTest(this,network,iterations)
            matrix = zeros(iterations+1,this.num_edge);
            matrix(1,:)=network;
            a = network<0.05;
            [mask,bigval] = this.my_largestCC(a);
            
%             this.clusterPlot(mask,this.all_edges(:, :,:),this.diagnosis);
            

            biggestArea = bigval;
            lcc=zeros(iterations,1);
            for i=1:iterations
                disp(i);
                inp2 = this.diagnosis;
                inp2 = inp2(randperm(length(inp2)));
                this.pval = zeros(this.num_edge,length(unique(inp2))-1);
                for j_edge = 1 : this.num_edge
                    inp1 = squeeze(this.all_edges(j_edge, :,:));
                    [d,p,~] = manova1(inp1,inp2);
                    this.dimensions(j_edge) = d;
                    this.pval(j_edge,:) = p(1);
                end
                matrix(i,:)=this.pval(:,1);
                a = this.pval(:,1)<0.05;
                [~,bigval] = this.my_largestCC(a);
                lcc(i)=bigval;
            end
            indices = find(sort(lcc)>=biggestArea);
%             disp(1-indices(1)/iterations);
        end
        function clusterPlot(this,mask,X,Y)
            [cx,cy] = find(mask);
            lindices = tril(ones(268));
            lindices(find(eye(268)))=0;
            lindices = find(lindices);
            P = zeros(length(unique(Y)),length(cx),length(X(1,1,:))); % 4*25*6
            %             X = squeeze(mean(X,2)); % take mean
            for g=1:length(unique(Y))
                indices = Y==g;
                Xg=squeeze(mean(X(:,indices,:),2));
                
                for i=1:length(cx)
                    for j=1:length( Xg(i,:))
                        network = zeros(268,268);
                        network(lindices) = Xg(:,j);
                        network = network'+network;
                        P(g,i,j)=network(cx(i),cy(i));
                    end
                end
            end
            %             subplot(1,3,3);
            d1=length(unique(Y));
            d2=5;%length(cx);
            for g=1:d1
                for j=1:d2
                    disp((g-1)*d2+j);
                    subplot(d1,d2,(g-1)*d2+j); %length(unique(Y))
                    bar(squeeze(P(g,j,:)))
                end
            end
            disp(P);
            
        end
    end
end

