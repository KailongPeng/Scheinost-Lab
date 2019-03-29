function summary = rCPM_pca(score,mapID)


x = load('/home/kailong/Scheinost-Lab/math/data/all_mats.mat');
% all_mats = x.all_mats(:,:,:,1);
% all_behav = sum_matrix_norm(:,1);
% all_behav = score(:,1);
% all_behav = common_factor;
% save('/home/kailong/Scheinost-Lab/math/data/all_behav.mat','all_behav')
p = 0.1;
field_rCPM_list = {'q_s', 'q_s_fold', 'r_pearson', 'r_rank', 'y', 'coef_total', 'coef0_total', 'lambda_total','p_pearson','p_rank'};
pc_list = [1:min([size(score,2),10])];
sum_performance= [];
for jj = 1:length(field_rCPM_list)
    field_rCPM = field_rCPM_list{jj};
    eval(['summary_' field_rCPM ' = [];'])
end
for curr_pc = 1:length(pc_list)
    pc = pc_list(curr_pc);
    all_behav = score(:,pc);
    for i = 1:2
        [q_s, q_s_fold, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total,p_pearson,p_rank] = ...
            siyuan_ridgeCPM(x.all_mats(:,:,mapID,:), all_behav,p);
        
        summary_q_s{curr_pc,i} = q_s;
        summary_q_s_fold{curr_pc,i} = q_s_fold;
        summary_r_pearson{curr_pc,i} = r_pearson;
        summary_r_rank{curr_pc,i} = r_rank;
        summary_y{curr_pc,i} = y;
        summary_coef_total{curr_pc,i} = coef_total;
        summary_coef0_total{curr_pc,i} = coef0_total;
        summary_lambda_total{curr_pc,i} = lambda_total;
        summary_p_pearson{curr_pc,i} = p_pearson;
        summary_p_rank{curr_pc,i} = p_rank;
        
%         for jj = 1:length(field_rCPM_list)
%             field_rCPM = field_rCPM_list{jj};
%             eval(['summary_' field_rCPM '{' num2str(i) ',' num2str(curr_pc) '} = ' field_rCPM ';'])
%         end
        %         sum_q_s{i,curr_pc} = q_s;
        
    end
    %     close all;
    %     figure;
    %     boxplot(tem)
    %     sum_performance{pc} = tem; tem = [];
    %     savefig([savedir 'performance of pc ' num2str(pc)])
end
for jj = 1:length(field_rCPM_list)
    field_rCPM = field_rCPM_list{jj};
    eval(['summary.' field_rCPM '= summary_' field_rCPM ';'])
end


% figure;
% boxplot(cell2mat(summary.q_s)');
% title(strrep('q_s PCA','_',' '))
% figure;
% boxplot(cell2mat(summary.r_pearson)');
% title(strrep('r_pearson PCA','_',' '))
% 
% figure;
% boxplot(cell2mat(summary.p_pearson)');
% title(strrep('p_pearson PCA','_',' '))
% 
% figure;
% boxplot(cell2mat(summary.r_rank)');
% title(strrep('r_rank PCA','_',' '))
% 
% figure;
% boxplot(cell2mat(summary.p_rank)');
% title(strrep('p_rank PCA','_',' '))
