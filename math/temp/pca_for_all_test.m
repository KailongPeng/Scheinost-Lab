clear all;
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
%check the correlation between different tests
[cr,~] = xcorr(remainTest_norm_no_nan_sum_matrix,0,'coeff');
cr = reshape(cr,[sqrt(length(cr)),sqrt(length(cr))]);
cr = cr - diag(diag(cr));
figure;imagesc(cr)

[noHighCor_table,noHighCor_table_VariableNames,noHighCor_norm_table] = excludeTooHighCorr(cr,remainTest_Variable_class,remainTest_VariableNames,remainTest_norm_no_nan_sum_matrix,remainTest_no_norm_no_nan_sum_all_test);
writetable(noHighCor_norm_table,'/home/kailong/Scheinost-Lab/math/data/noHighCor_norm_SelectedTest_test_no_norm_no_nan_header','Delimiter',',')
writetable(noHighCor_table,'/home/kailong/Scheinost-Lab/math/data/noHighCor_SelectedTest_test_no_norm_no_nan_header','Delimiter',',')
writetable(noHighCor_table_VariableNames,'/home/kailong/Scheinost-Lab/math/data/SelectedTest_VariableNames','Delimiter',',')

%%
%pca
% [pca_coefficient,score,latent] = pca(norm_no_nan_sum_math_matrix,'algorithm','als');
[pca_coefficient,score,latent] = pca(norm_no_nan_sum_all_test);
%for pruning PC
varianceThreshold = 0.5;
minComponent = 1;
while sum(latent(1:minComponent))/sum(latent(:)) < varianceThreshold %threshold for variance to determine minimum amount of principal components
    minComponent = minComponent+1;
end

summary = rCPM_pca(score,all_mapID);
savedir = '/home/kailong/Scheinost-Lab/math/plot/rCPM/pca_no_nan/pearson/pca/all_test/';
if ~isdir(savedir); mkdir(savedir); end
save([savedir 'performance of pc'],'summary');
1;
