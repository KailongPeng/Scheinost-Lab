% use Siyuan's new code to analyse the latent factor scores

clear all;close all;clc;
% path = '/home/kailong/Scheinost-Lab/math/data/LatentFactorEstimate/';
path = '/home/kailong/Scheinost-Lab/math/data/Norm_LatentFactorEstimate/';
LatentFactorList = dir([path '*.mat']);
LatentFactorList = kailong_extractfield(LatentFactorList,'name');
all_mats = load('/home/kailong/Scheinost-Lab/math/data/all_mats.mat');
all_mats = all_mats.all_mats;
load('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header','all_mapID')
all_mats = all_mats(:,:,all_mapID,:);

q_s = [];
r_pearson = [];
r_rank = [];
y = [];
new_behav = [];
all_edge_weight = [];
all_behav_weight = [];
all_task_weight = [];
lambda_total = [];
for curr_LatentFactor = 3:size(LatentFactorList,2)
    fprintf('only doing 3 factors\n');pause;
    LatentFactor = LatentFactorList{curr_LatentFactor};
    all_behav = [];
    all_behav = load([path LatentFactor]);
    all_behav = all_behav.data;
    thresh1 = 0.3;
    thresh2 = 0.1;
    tStart1 = tic;
    lambda =[];
    k = 10;
    [q_s{curr_LatentFactor}, r_pearson{curr_LatentFactor}, r_rank{curr_LatentFactor}, y{curr_LatentFactor}, new_behav{curr_LatentFactor}, all_edge_weight{curr_LatentFactor}, all_behav_weight{curr_LatentFactor}, all_task_weight{curr_LatentFactor}, lambda_total{curr_LatentFactor}] = ...
        mmCPM_ridge(all_mats, all_behav, thresh1, thresh2, lambda, k);
    tElapsed = toc(tStart1)
    temp = [];
    temp = q_s{curr_LatentFactor};
    temp(temp<0) = 0;
    temp = sqrt(temp);
    q_s{curr_LatentFactor} = temp;
    1;
end

AllQs = cell2mat(q_s);
AllQs(AllQs<0) = nan;
AllQs = sqrt(AllQs);
figure;boxplot(AllQs);xlabel('latent factor number');ylabel('sqrt(qs)');
AllRPearson = cell2mat(r_pearson);
figure;boxplot(AllRPearson);xlabel('latent factor number');ylabel('r pearson');
AllRRank = cell2mat(r_rank);
figure;boxplot(AllRRank);xlabel('latent factor number');ylabel('r rank');


AllY = cell2mat(y);
figure;boxplot(AllY);xlabel('latent factor number');ylabel('y');

figure;bar3(AllQs);xlabel('latent factor number');ylabel('kfold');zlabel('sqrt(qs)');
figure;bar3(AllRPearson);xlabel('latent factor number');ylabel('kfold');zlabel('r pearson');
figure;bar3(AllRRank);xlabel('latent factor number');ylabel('kfold');zlabel('r rank');
%%
clear all;close all;clc;
addpath(genpath('/home/kailong/Scheinost-Lab'))
% path = '/home/kailong/Scheinost-Lab/math/data/Norm_LatentFactorEstimate/';
% LatentFactorList = dir([path '*.mat']);
% LatentFactorList = kailong_extractfield(LatentFactorList,'name');
all_mats = load('/home/kailong/Scheinost-Lab/math/data/all_mats.mat');
all_mats = all_mats.all_mats;
load('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header','all_mapID')
all_mats = all_mats(:,:,all_mapID,:);

load('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header')
no_norm_no_nan_sum_all_test = table2array(table);
norm_no_nan_sum_all_test = normalize(table2array(table),1);
%pick only the test selected here 
% allRemainTestList = ["ADHD_Total_%","AWMA-S_VerbalSTM_StS","AWMA-S_VerbalWM_StS","AWMA-S_VisuoSpatialSTM_StS","AWMA-S_VisuoSpatialWM_StS","CMAT_BasicCalc_Comp_Quotient","CTOPP_EL_StS","CTOPP_BW_StS","CTOPP_PhonAwareness_Comp","CTOPP_RapidNaming_Comp","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA-2_Attitudes_StS","TOWRE_Total_StS","WASI_Vocab_T-Score","WASI_BD_T-Score","WASI_Sim_T-Score","WASI_MR_T-Score","WASI_VIQ","WASI_PIQ","WASI_FSIQ","WJ-III_WordID_StS","WJ-III_WA_StS","WJ-III_PassComp_StS","WJ-III_MathFluency_StS","WJ-III_SpatialRelations_StS","WJ-III_BRS_StS"];
allRemainTestList = ["ADHD_Total_0x25","AWMA0x2DS_VerbalSTM_StS","AWMA0x2DS_VerbalWM_StS","AWMA0x2DS_VisuoSpatialSTM_StS","AWMA0x2DS_VisuoSpatialWM_StS","CMAT_BasicCalc_Comp_Quotient","CTOPP_EL_StS","CTOPP_BW_StS","CTOPP_PhonAwareness_Comp","CTOPP_RapidNaming_Comp","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA0x2D2_Attitudes_StS","TOWRE_Total_StS","WASI_Vocab_T0x2DScore","WASI_BD_T0x2DScore","WASI_Sim_T0x2DScore","WASI_MR_T0x2DScore","WASI_VIQ","WASI_PIQ","WASI_FSIQ","WJ0x2DIII_WordID_StS","WJ0x2DIII_WA_StS","WJ0x2DIII_PassComp_StS","WJ0x2DIII_MathFluency_StS","WJ0x2DIII_SpatialRelations_StS","WJ0x2DIII_BRS_StS"];
remainTestID = cellfun(@(x) sum(contains(allRemainTestList,x)),char_VariableNames);
remainTestID = find(remainTestID==1);
char_VariableNames = char(char_VariableNames);
remainTest_VariableNames = char_VariableNames(remainTestID,:);
remainTest_VariableNames = cellstr(remainTest_VariableNames);
remainTest_no_norm_no_nan_sum_all_test = no_norm_no_nan_sum_all_test(:,remainTestID);
remainTest_norm_no_nan_sum_matrix = normalize(remainTest_no_norm_no_nan_sum_all_test,1);
remainTest_Variable_class = Variable_class(remainTestID);

all_behav = remainTest_no_norm_no_nan_sum_all_test;

save('/home/kailong/Scheinost-Lab/code_from_javid/rCPM/input','all_behav','all_mats')

thresh1 = 0.3;
thresh2 = 0.1;
tStart1 = tic;
lambda =[];
k = 20;

LatentFactorList = [1:5];%size(all_behav,2)];
for curr_LatentFactor = 1:size(LatentFactorList,2)
%     [q_s, r_pearson, r_rank, y, new_behav, all_edge_weight, all_behav_weight, all_task_weight, lambda_total] = ...
%         kailong_mmCPM_ridge(all_mats, all_behav, thresh1, thresh2, lambda, k);
    numOfFactor = curr_LatentFactor;
%     numOfPC = 2;
    [q_s{curr_LatentFactor}, r_pearson{curr_LatentFactor}, r_rank{curr_LatentFactor}, ...
        y{curr_LatentFactor}, new_behav{curr_LatentFactor}, all_edge_weight{curr_LatentFactor}, ...
        all_behav_weight{curr_LatentFactor}, all_task_weight{curr_LatentFactor}, ...
        lambda_total{curr_LatentFactor},FA_Lambda{curr_LatentFactor}] = ...
        kailong_mmCPM_ridge(all_mats, all_behav, thresh1, thresh2, lambda, k, numOfFactor, numOfPC);
end
tElapsed = toc(tStart1)
save('/home/kailong/Scheinost-Lab/math/FactorAnalysis/workspace')

figure;
for curr_LatentFactor = 1:size(LatentFactorList,2)
    clear temp
    temp = FA_Lambda{curr_LatentFactor}{1};
    for i_fold = 2:k
        temp = temp + FA_Lambda{curr_LatentFactor}{i_fold};
    end
    ave_FA_Lambda{curr_LatentFactor} = temp/k;
    
    FA_Lambda_mat{curr_LatentFactor} = (ave_FA_Lambda{curr_LatentFactor});
    FA_Lambda_mat{curr_LatentFactor}(FA_Lambda_mat{curr_LatentFactor}<0.5) = nan;
    subplot(2,3,curr_LatentFactor)
    imagesc(FA_Lambda_mat{curr_LatentFactor});
    yticks([1:size(allRemainTestList,2)])
    yticklabels(strrep(allRemainTestList,'_','-'))
end

AllQs = cell2mat(q_s);
AllQs(AllQs<0) = nan;
AllQs = sqrt(AllQs);
figure;subplot(2,3,1)
boxplot(AllQs);xlabel('latent factor number');ylabel('sqrt(qs)');
AllRPearson = cell2mat(r_pearson);
subplot(2,3,2)
boxplot(AllRPearson);xlabel('latent factor number');ylabel('r pearson');
AllRRank = cell2mat(r_rank);
subplot(2,3,3)
boxplot(AllRRank);xlabel('latent factor number');ylabel('r rank');

subplot(2,3,4)
bar3(AllQs);xlabel('latent factor number');ylabel('kfold');zlabel('sqrt(qs)');
subplot(2,3,5)
bar3(AllRPearson);xlabel('latent factor number');ylabel('kfold');zlabel('r pearson');
subplot(2,3,6)
bar3(AllRRank);xlabel('latent factor number');ylabel('kfold');zlabel('r rank');

AllY = cell2mat(y);
figure;boxplot(AllY);xlabel('latent factor number');ylabel('y');

