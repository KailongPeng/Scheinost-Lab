% try some positive control because latent factor analysis has failed.
% try some standard metrics like vocabulary, matrix reasoning, verbal IQ,
% performance IQ, and full IQ

% function PredictIntelligence()

clear all;close all;clc;

%tests that I want to try:
% IntelligenceTestList = ["WASI_Vocab_T0x2DScore","WASI_BD_T0x2DScore","WASI_Sim_T0x2DScore","WASI_MR_T0x2DScore","WASI_VIQ","WASI_PIQ","WASI_FSIQ"];
IntelligenceTestList = ["ADHD_Total_0x25","AWMA0x2DS_VerbalSTM_StS","AWMA0x2DS_VerbalWM_StS","AWMA0x2DS_VisuoSpatialSTM_StS","AWMA0x2DS_VisuoSpatialWM_StS","CMAT_BasicCalc_Comp_Quotient","CTOPP_EL_StS","CTOPP_BW_StS","CTOPP_PhonAwareness_Comp","CTOPP_RapidNaming_Comp","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA0x2D2_Attitudes_StS","TOWRE_Total_StS","WASI_Vocab_T0x2DScore","WASI_BD_T0x2DScore","WASI_Sim_T0x2DScore","WASI_MR_T0x2DScore","WASI_VIQ","WASI_PIQ","WASI_FSIQ","WJ0x2DIII_WordID_StS","WJ0x2DIII_WA_StS","WJ0x2DIII_PassComp_StS","WJ0x2DIII_MathFluency_StS","WJ0x2DIII_SpatialRelations_StS","WJ0x2DIII_BRS_StS"];

% remainTestList = ["CMAT_BasicCalc_Comp_Quotient","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA0x2D2_Attitudes_StS"];
% remainTestID = cellfun(@(x) sum(contains(remainTestList,x)),char_math_VariableNames);
% remainTestID = find(remainTestID==1);
output = [];
complete_all_mats = load('/home/kailong/Scheinost-Lab/math/data/all_mats.mat');
All_bad_id = load('/home/kailong/Scheinost-Lab/math/data/motion.mat', 'motion_para', 'bad_id');
All_data_all_test = load('/home/kailong/Scheinost-Lab/math/data/PredictIntelligenceTests/data_all_test','VariableNames',...
    'data_all_test');
% parfor (curr_IntelligenceTest = 1:size(IntelligenceTestList,2),4)
for curr_IntelligenceTest = 1:size(IntelligenceTestList,2)
    IntelligenceTest = IntelligenceTestList(curr_IntelligenceTest);
    VariableNames = All_data_all_test.VariableNames;
    char_VariableNames = char(VariableNames);
    
    remainTestID = find(cellfun(@(x) sum(contains(IntelligenceTest,x)),VariableNames)==1);
    
    data_all_test = All_data_all_test.data_all_test;
    all_mapID = [1:size(data_all_test,1)];
    
    %deal with bad subjects -high motion
    bad_id = All_bad_id.bad_id;
    data_all_test(bad_id,:) = [];%data_all_test original 132subjects*120tests %now 84sub*120tests
    all_mapID(bad_id) = [];
    
    data_all_test = data_all_test(:,remainTestID); %now 84sub*1tests
    all_mapID = all_mapID(~isnan(data_all_test));
    data_all_test = data_all_test(~isnan(data_all_test));
    all_behav = data_all_test;
    all_mats = complete_all_mats.all_mats; % original: 268node*268node*132sub*7task
    all_mats = all_mats(:,:,all_mapID,:); % new: 268node*268node*84*7task
    
    thresh = 0.1;
    v_alpha = [];
    lambda = [];
    k = 5;
    [q_s, q_s_fold, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total,p_pearson,p_rank] = ...
        siyuan_ridgeCPM(all_mats, all_behav, thresh, v_alpha, lambda, k);
    output{curr_IntelligenceTest}.q_s = q_s;
    output{curr_IntelligenceTest}.q_s_fold = q_s_fold;
    output{curr_IntelligenceTest}.r_pearson = r_pearson;
    output{curr_IntelligenceTest}.r_rank = r_rank;
    output{curr_IntelligenceTest}.y = y;
    output{curr_IntelligenceTest}.coef_total = coef_total;
    output{curr_IntelligenceTest}.coef0_total = coef0_total;
    output{curr_IntelligenceTest}.lambda_total = lambda_total;
    output{curr_IntelligenceTest}.p_pearson = p_pearson;
    output{curr_IntelligenceTest}.p_rank = p_rank;
end
save(['/home/kailong/Scheinost-Lab/math/temp/output-' datestr(datetime('now'))],'output')

figure;
q_s = cell2mat(kailong_extractfield(output,'q_s'));
subplot(5,1,1)
plot(q_s,'or')
xlabel('q s')
set(gca,'XTick',[])
q = q_s;q(q<0) = nan;q = sqrt(q);
r_pearson = cell2mat(kailong_extractfield(output,'r_pearson'));
subplot(5,1,2)
plot(r_pearson,'or')
xlabel('r pearson')
set(gca,'XTick',[])
r_rank = cell2mat(kailong_extractfield(output,'r_rank'));
subplot(5,1,3)
plot(r_rank,'or')
xlabel('r rank')
set(gca,'XTick',[])
p_pearson = cell2mat(kailong_extractfield(output,'p_pearson'));
subplot(5,1,4)
plot(p_pearson,'or')
xlabel('p pearson')
set(gca,'XTick',[])
p_rank = cell2mat(kailong_extractfield(output,'p_rank'));
subplot(5,1,5)
plot(p_rank,'or')
xlabel('p rank')
set(gca,'XTick',[])
xticks([1:size(IntelligenceTestList,2)])
xticklabels(strrep(IntelligenceTestList,'_','-'));
xtickangle(90)

temp = r_pearson;
temp(temp<0.1) = nan;
RefinedIntelligenceTestList = IntelligenceTestList(1,(~isnan(temp) == 1));
save('RefinedIntelligenceTestList01','RefinedIntelligenceTestList')
% end


%%
% function UseHighPrectingPowerPairToDoFactorAnalysisMmCPM()
clear all;close all;clc;
load('/home/kailong/Scheinost-Lab/math/temp/RefinedIntelligenceTestList01')
output = [];
complete_all_mats = load('/home/kailong/Scheinost-Lab/math/data/all_mats.mat');
All_bad_id = load('/home/kailong/Scheinost-Lab/math/data/motion.mat', 'motion_para', 'bad_id');
All_data_all_test = load('/home/kailong/Scheinost-Lab/math/data/PredictIntelligenceTests/data_all_test','VariableNames',...
    'data_all_test');
% parfor (curr_IntelligenceTest = 1:size(IntelligenceTestList,2),4)
% for curr_IntelligenceTest = 1:size(IntelligenceTestList,2)
VariableNames = All_data_all_test.VariableNames;
char_VariableNames = char(VariableNames);

remainTestID = (cellfun(@(x) sum(contains(RefinedIntelligenceTestList,x)),VariableNames)==1);

data_all_test = All_data_all_test.data_all_test;
all_mapID = [1:size(data_all_test,1)];

%deal with bad subjects -high motion
bad_id = All_bad_id.bad_id;
data_all_test(bad_id,:) = [];%data_all_test original 132subjects*120tests %now 84sub*120tests
all_mapID(bad_id) = [];

data_all_test = data_all_test(:,remainTestID); %now 84sub*7tests
all_mapID = all_mapID(~sum(isnan(data_all_test),2));%now 78sub*7tests
data_all_test = data_all_test(~sum(isnan(data_all_test),2),:);%now 78sub*7tests
all_behav = data_all_test;
all_mats = complete_all_mats.all_mats; % original: 268node*268node*132sub*7task
all_mats = all_mats(:,:,all_mapID,:); % new: 268node*268node*78*7task
%
% thresh = 0.1;
% v_alpha = [];
% lambda = [];
% k = 5;
thresh1 = 0.3;
thresh2 = 0.1;
tStart1 = tic;
lambda =[];
k = 5;
for curr_LatentFactor = 1:3%1:size(LatentFactorList,2)
    numOfFactor = curr_LatentFactor;
    numOfPC = [];
    singleFactor = [];
    %     [q_s, q_s_fold, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total,p_pearson,p_rank] = ...
    %         siyuan_ridgeCPM(all_mats, all_behav, thresh, v_alpha, lambda, k);
    [q_s, r_pearson, r_rank, y, new_behav, all_edge_weight, all_behav_weight, all_task_weight, lambda_total,FA_Lambda,PCA_Lambda] = ...
        kailong_mmCPM_ridge(all_mats, all_behav, thresh1, thresh2, lambda, k, numOfFactor, numOfPC,singleFactor);
    
    output{curr_LatentFactor}.q_s = q_s;
    %     output{curr_LatentFactor}.q_s_fold = q_s_fold;
    output{curr_LatentFactor}.r_pearson = r_pearson;
    output{curr_LatentFactor}.r_rank = r_rank;
    output{curr_LatentFactor}.y = y;
    output{curr_LatentFactor}.new_behav = new_behav;
    output{curr_LatentFactor}.all_edge_weight = all_edge_weight;
    output{curr_LatentFactor}.all_behav_weight = all_behav_weight;
    output{curr_LatentFactor}.all_task_weight = all_task_weight;
    output{curr_LatentFactor}.lambda_total = lambda_total;
    output{curr_LatentFactor}.FA_Lambda = FA_Lambda;
    output{curr_LatentFactor}.PCA_Lambda = PCA_Lambda;
    %     output{curr_LatentFactor}.coef_total = coef_total;
    %     output{curr_LatentFactor}.coef0_total = coef0_total;
    %     output{curr_LatentFactor}.p_pearson = p_pearson;
    %     output{curr_LatentFactor}.p_rank = p_rank;
end
save(['/home/kailong/Scheinost-Lab/math/temp/output-' datestr(datetime('now'))],'output')


figure;
q_s = mean(cell2mat(kailong_extractfield(output,'q_s')));
subplot(4,1,1)
plot(q_s,'or')
xlabel('q s')
set(gca,'XTick',[])
q = cell2mat(kailong_extractfield(output,'q_s'));q(q<0) = nan;q = nanmean(sqrt(q));
subplot(4,1,2)
plot(q,'or')
xlabel('q')
set(gca,'XTick',[])
r_pearson = mean(cell2mat(kailong_extractfield(output,'r_pearson')));
subplot(4,1,3)
plot(r_pearson,'or')
xlabel('r pearson')
set(gca,'XTick',[])
r_rank = mean(cell2mat(kailong_extractfield(output,'r_rank')));
subplot(4,1,4)
plot(r_rank,'or')
xlabel('r rank')
set(gca,'XTick',[])
x_label = ["1 latent factor" "2 latent factors" "3 latent factors"];
xticks([1:size(x_label,2)])
xticklabels(strrep(x_label,'_','-'));
xtickangle(45)

%%
% load('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header','all_mapID')
% all_mats = all_mats(:,:,all_mapID,:);
% load('/home/kailong/Scheinost-Lab/math/data/motion.mat', 'motion_para', 'bad_id');
%
% load('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header')
% no_norm_no_nan_sum_all_test = table2array(table);
% norm_no_nan_sum_all_test = normalize(table2array(table),1);
% %pick only the test selected here
% % allRemainTestList = ["ADHD_Total_%","AWMA-S_VerbalSTM_StS","AWMA-S_VerbalWM_StS","AWMA-S_VisuoSpatialSTM_StS","AWMA-S_VisuoSpatialWM_StS","CMAT_BasicCalc_Comp_Quotient","CTOPP_EL_StS","CTOPP_BW_StS","CTOPP_PhonAwareness_Comp","CTOPP_RapidNaming_Comp","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA-2_Attitudes_StS","TOWRE_Total_StS","WASI_Vocab_T-Score","WASI_BD_T-Score","WASI_Sim_T-Score","WASI_MR_T-Score","WASI_VIQ","WASI_PIQ","WASI_FSIQ","WJ-III_WordID_StS","WJ-III_WA_StS","WJ-III_PassComp_StS","WJ-III_MathFluency_StS","WJ-III_SpatialRelations_StS","WJ-III_BRS_StS"];
% allRemainTestList = ["ADHD_Total_0x25","AWMA0x2DS_VerbalSTM_StS","AWMA0x2DS_VerbalWM_StS","AWMA0x2DS_VisuoSpatialSTM_StS","AWMA0x2DS_VisuoSpatialWM_StS","CMAT_BasicCalc_Comp_Quotient","CTOPP_EL_StS","CTOPP_BW_StS","CTOPP_PhonAwareness_Comp","CTOPP_RapidNaming_Comp","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA0x2D2_Attitudes_StS","TOWRE_Total_StS","WASI_Vocab_T0x2DScore","WASI_BD_T0x2DScore","WASI_Sim_T0x2DScore","WASI_MR_T0x2DScore","WASI_VIQ","WASI_PIQ","WASI_FSIQ","WJ0x2DIII_WordID_StS","WJ0x2DIII_WA_StS","WJ0x2DIII_PassComp_StS","WJ0x2DIII_MathFluency_StS","WJ0x2DIII_SpatialRelations_StS","WJ0x2DIII_BRS_StS"];
% remainTestID = cellfun(@(x) sum(contains(allRemainTestList,x)),char_VariableNames);
% remainTestID = find(remainTestID==1);
% char_VariableNames = char(char_VariableNames);
% remainTest_VariableNames = char_VariableNames(remainTestID,:);
% remainTest_VariableNames = cellstr(remainTest_VariableNames);
% remainTest_no_norm_no_nan_sum_all_test = no_norm_no_nan_sum_all_test(:,remainTestID);
% remainTest_norm_no_nan_sum_matrix = normalize(remainTest_no_norm_no_nan_sum_all_test,1);
% remainTest_Variable_class = Variable_class(remainTestID);
%
% all_behav = remainTest_no_norm_no_nan_sum_all_test;
