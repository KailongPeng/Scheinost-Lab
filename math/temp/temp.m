restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab'));

load('/home/kailong/Scheinost-Lab/math/data/math_test')
sum_matrix_norm = normalize(sum_matrix,1);
mapID = [1:size(sum_matrix_norm,1)];
mapID = mapID(all(~isnan(sum_matrix_norm),2));
sum_matrix_norm_no_nan = sum_matrix_norm(all(~isnan(sum_matrix_norm),2),:);
save('/home/kailong/Scheinost-Lab/math/data/math_test_norm_no_nan','sum_matrix_norm_no_nan')
csvwrite('/home/kailong/Scheinost-Lab/math/data/math_test_norm_no_nan',sum_matrix_norm_no_nan)
%check the correlation between different tests
[cr,~] = xcorr(sum_matrix_norm_no_nan,0,'coeff');
cr = reshape(cr,[17,17]);
cr = cr - diag(diag(cr));
figure
imagesc(cr)

%pca
% [pca_coefficient,score,latent] = pca(sum_matrix_norm_no_nan,'algorithm','als');
[pca_coefficient,score,latent] = pca(sum_matrix_norm_no_nan);
%for pruning PC
varianceThreshold = 0.5;
minComponent = 1;
while sum(latent(1:minComponent))/sum(latent(:)) < varianceThreshold %threshold for variance to determine minimum amount of principal components
    minComponent = minComponent+1;
end

% use_scores(:,:) = score(:,1:minComponent);
% nanmean(use_scores);
% nanvar(use_scores);
% pca_coefficient(:,1);

%factor analysis

for num_common_factor = 1:10
    [lambda,psi,T,stats,F] = factoran(sum_matrix_norm_no_nan,num_common_factor);
    if stats.p < 0.05
        break;
    end
end
common_factor = repmat(lambda',[size(sum_matrix_norm_no_nan,1),1]).*sum_matrix_norm_no_nan;
common_factor = sum(common_factor,2);
% figure;plot([1:17],lambda,'r.');hold on;plot([1:17],1-psi,'b.');plot([1:17],F,'g.');


x = load('/home/kailong/Scheinost-Lab/math/data/all_mats.mat');
all_mats = x.all_mats(:,:,:,1);
% all_behav = sum_matrix_norm(:,1);
% all_behav = score(:,1);
% all_behav = common_factor;
% save('/home/kailong/Scheinost-Lab/math/data/all_behav.mat','all_behav')
p = 0.1;
kfold = 10;
% corr_type = 'spearman';
corr_type = 'pearson';
LinearFlag = 1;

clear sum_performance
savedir = '/home/kailong/Scheinost-Lab/math/plot/rCPM/pca_no_nan/pearson/pca/';
if ~isdir(savedir); mkdir(savedir); end
pc_list = [1 2 10 17];
for curr_pc = 1:length(pc_list)
    pc = pc_list(curr_pc);
    all_behav = score(:,pc);
    tem = [];
    for i = 1:10
        [q_s, q_s_fold, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total] = ...
            siyuan_ridgeCPM(x.all_mats(:,:,mapID,:), all_behav,p);
%         [y_predict,performance] = cpm_main(all_mats(:,:,mapID),all_behav,corr_type,LinearFlag,'pthresh',p,'kfolds',kfold);
        field_rCPM_list = {'q_s', 'q_s_fold', 'r_pearson', 'r_rank', 'y', 'coef_total', 'coef0_total', 'lambda_tota'};
        for jj = 1:length(field_rCPM_list)
            field_rCPM = field_rCPM_list{jj};
            eval(['summary_' field_rCPM '{' i ',' curr_pc '} = ' field_rCPM ';'])
        end
%         sum_q_s{i,curr_pc} = q_s;
        
    end
%     close all;
%     figure;
%     boxplot(tem)
%     sum_performance{pc} = tem; tem = [];
%     savefig([savedir 'performance of pc ' num2str(pc)])
end
save([savedir 'performance of pc'],'summary')
figure;
boxplot(cell2mat(summary.r_rank))