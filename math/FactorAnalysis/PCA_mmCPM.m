% use Siyuan's new code to analyse the PCA scores

clear all;close all;clc;
% path = '/home/kailong/Scheinost-Lab/math/data/LatentFactorEstimate/';
% path = '/home/kailong/Scheinost-Lab/math/data/Norm_LatentFactorEstimate/';
% LatentFactorList = dir([path '*.mat']);
% LatentFactorList = kailong_extractfield(LatentFactorList,'name');
all_mats = load('/home/kailong/Scheinost-Lab/math/data/all_mats.mat');
all_mats = all_mats.all_mats;
load('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header')%,'all_mapID','table')
all_mats = all_mats(:,:,all_mapID,:);

no_norm_no_nan_sum_all_test = table2array(table);
norm_no_nan_sum_all_test = normalize(table2array(table),1);

allRemainTestList = ["ADHD_Total_0x25","AWMA0x2DS_VerbalSTM_StS","AWMA0x2DS_VerbalWM_StS","AWMA0x2DS_VisuoSpatialSTM_StS","AWMA0x2DS_VisuoSpatialWM_StS","CMAT_BasicCalc_Comp_Quotient","CTOPP_EL_StS","CTOPP_BW_StS","CTOPP_PhonAwareness_Comp","CTOPP_RapidNaming_Comp","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA0x2D2_Attitudes_StS","TOWRE_Total_StS","WASI_Vocab_T0x2DScore","WASI_BD_T0x2DScore","WASI_Sim_T0x2DScore","WASI_MR_T0x2DScore","WASI_VIQ","WASI_PIQ","WASI_FSIQ","WJ0x2DIII_WordID_StS","WJ0x2DIII_WA_StS","WJ0x2DIII_PassComp_StS","WJ0x2DIII_MathFluency_StS","WJ0x2DIII_SpatialRelations_StS","WJ0x2DIII_BRS_StS"];
remainTestID = cellfun(@(x) sum(contains(allRemainTestList,x)),char_VariableNames);
remainTestID = find(remainTestID==1);
char_VariableNames = char(char_VariableNames);
remainTest_VariableNames = char_VariableNames(remainTestID,:);
remainTest_VariableNames = cellstr(remainTest_VariableNames);
remainTest_no_norm_no_nan_sum_all_test = no_norm_no_nan_sum_all_test(:,remainTestID);
% remainTest_norm_no_nan_sum_matrix = normalize(remainTest_no_norm_no_nan_sum_all_test,1);

% [pca_coefficient,score,latent] = pca(remainTest_norm_no_nan_sum_matrix,'algorithm','als');

q_s = [];
r_pearson = [];
r_rank = [];
y = [];
new_behav = [];
all_edge_weight = [];
all_behav_weight = [];
all_task_weight = [];
lambda_total = [];
for curr_pc = 1:size(score,2)
    all_behav = [];
%     all_behav = score(:,curr_pc);
    all_behav = remainTest_no_norm_no_nan_sum_all_test;
    thresh1 = 0.3;
    thresh2 = 0.1;
    tStart1 = tic;
    lambda =[];
    k = 10;
    [q_s{curr_pc}, r_pearson{curr_pc}, r_rank{curr_pc}, y{curr_pc}, new_behav{curr_pc}, all_edge_weight{curr_pc}, all_behav_weight{curr_pc}, all_task_weight{curr_pc}, lambda_total{curr_pc}] = ...
        kailong_mmCPM_ridge(all_mats, all_behav, thresh1, thresh2, lambda, k);
    tElapsed = toc(tStart1)
    temp = [];
    temp = q_s{curr_pc};
    temp(temp<0) = 0;
    temp = sqrt(temp);
    q_s{curr_pc} = temp;
    1;
end

AllQs = cell2mat(q_s);
AllQs(AllQs==0) = nan;
figure;boxplot(AllQs);xlabel('PC number');ylabel('sqrt(qs)');
AllRPearson = cell2mat(r_pearson);
figure;boxplot(AllRPearson);xlabel('PC number');ylabel('r pearson');
AllRRank = cell2mat(r_rank);
figure;boxplot(AllRRank);xlabel('PC number');ylabel('r rank');


AllY = cell2mat(y);
figure;boxplot(AllY);xlabel('PC number');ylabel('y');

figure;bar3(AllQs);xlabel('PC number');ylabel('kfold');zlabel('sqrt(qs)');
figure;bar3(AllRPearson);xlabel('PC number');ylabel('kfold');zlabel('r pearson');
figure;bar3(AllRRank);xlabel('PC number');ylabel('kfold');zlabel('r rank');
1;