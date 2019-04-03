field_rCPM_list = {'q_s', 'q_s_fold', 'r_pearson', 'r_rank', 'y', 'coef_total', 'coef0_total', 'lambda_total','p_pearson','p_rank'};

pc_list = [1:size(score,2)];
for curr_pc = 1:length(pc_list)
    pc = pc_list(curr_pc);
    all_behav = score(:,pc);
    tem = [];
    for i = 1:2
%         [q_s, q_s_fold, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total,p_pearson,p_rank] = ...
%             siyuan_ridgeCPM(x.all_mats(:,:,mapID,:), all_behav,p);
        for jj = 1:length(field_rCPM_list)
            field_rCPM = field_rCPM_list{jj};
            fprintf(['summary_' field_rCPM '{' num2str(i) ',' num2str(curr_pc) '} = ' field_rCPM ';\n'])
        end
        %         sum_q_s{i,curr_pc} = q_s;
        
    end
    %     close all;
    %     figure;
    %     boxplot(tem)
    %     sum_performance{pc} = tem; tem = [];
    %     savefig([savedir 'performance of pc ' num2str(pc)])
end