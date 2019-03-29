clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab'));
cd('/home/kailong/Scheinost-Lab')

% %dustin data： functional and anatomical data
% %check for motion and the registration
% cd('/data_dustin/math2/results/')
% load('/data_dustin/math2/results/');
% 
% %siyuan data
% cd('/mnt/store4/mri_group/siyuan_data/math');
% load('/mnt/store4/mri_group/siyuan_data/math');

%%
%behavioral measures
path = '/data_dustin/math2/phenotype/';
cd(path);
fnameList = dir([path '*.json']);
fnameList = kailong_extractfield(fnameList,'name');
selected_task_list = {'cmat' 'keymath-3' 'toma-2'};
clear task_description
for curr_task = 1:length(fnameList)
    fname = fnameList{curr_task};
    taskname = kailong_minus(fname,'.json');

    filename = ['ses-T1/' taskname '.tsv'];
    try task_data(curr_task).sesT1 = tdfread(filename,'tab');end
    filename = ['ses-T2/' taskname '.tsv'];
    try task_data(curr_task).sesT2 = tdfread(filename); end
    task_data(curr_task).taskname = taskname;%'adhd-rs.json';
    curr_math_index = find(cellfun((@(x) strcmp(x,taskname)),selected_task_list)==1);
    if ~isempty(curr_math_index)
%         curr_math_index = find(cellfun((@(x) contains(x,taskname)),selected_task_list)==1);
        task_math(curr_math_index).sesT1 = task_data(curr_task).sesT1;
        task_math(curr_math_index).sesT2 = task_data(curr_task).sesT2;
        task_math(curr_math_index).taskname = taskname;%'adhd-rs.json';
    end
    
    task_description(curr_task).task = jsondecode(fileread(fname));%'adhd-rs.json';
    task_description(curr_task).taskname = taskname;%'adhd-rs.json';
    filename = [];
    fname = [];
    taskname = [];
end

% for curr_task = 1:size(task_description,2)
%     VariableNames = fields(task_description(1).task)
% end

% task_math(1).sesT1.participant_id_num = str2num(task_math(1).sesT1.participant_id(end-3,end));
% temp = task_math(3).sesT1.TOMA0x2D2_Attitudes_StS;
% temp2 = get_number(temp);
sum_math_matrix = [];
math_VariableNames = [];
curr_VariableNames = 1;
curr_class = 1;
for curr_task = 1:length(task_math)
    curr_task_fields_list = fields(task_math(curr_task).sesT1);
    for curr_task_field = 1:length(curr_task_fields_list)
        task_field = curr_task_fields_list{curr_task_field};
        if ~contains(task_field,'participant_id')
            if strcmp(class(eval(['task_math(curr_task).sesT1.' task_field ''])),'char')
                temp = [];
                temp = eval(['get_number(task_math(curr_task).sesT1.' task_field ')']);
                sum_math_matrix = [sum_math_matrix temp(:)];
            else
                temp = [];
                temp = eval(['task_math(curr_task).sesT1.' task_field ';']);
                sum_math_matrix = [sum_math_matrix temp(:)];
            end
            math_Variable_class(curr_VariableNames) = curr_class;
            math_VariableNames{curr_VariableNames} = task_field;
            curr_VariableNames = curr_VariableNames + 1;
        end
    end
    curr_class = curr_class + 1;
end
char_math_VariableNames = char(math_VariableNames);
char_math_VariableNames = char_math_VariableNames(logical(nansum(sum_math_matrix)),:);
char_math_VariableNames = cellstr(char_math_VariableNames);
no_norm_no_nan_sum_math_matrix = sum_math_matrix(:,logical(nansum(sum_math_matrix)));

math_mapID = [1:size(no_norm_no_nan_sum_math_matrix,1)];
math_mapID = math_mapID(all(~isnan(no_norm_no_nan_sum_math_matrix),2));
no_norm_no_nan_sum_math_matrix = no_norm_no_nan_sum_math_matrix(all(~isnan(no_norm_no_nan_sum_math_matrix),2),:);
math_table = [];
math_table = array2table(no_norm_no_nan_sum_math_matrix,'VariableNames',char_math_VariableNames);
% writetable(math_table,'/home/kailong/Scheinost-Lab/math/data/math_test_no_norm_no_nan_header','Delimiter',',')
save('/home/kailong/Scheinost-Lab/math/data/math_test_no_norm_no_nan_header','math_table','math_mapID','char_math_VariableNames','math_Variable_class')

table_math_VariableNames = cell2table(char_math_VariableNames');
% writetable(table_math_VariableNames,'/home/kailong/Scheinost-Lab/math/data/math_VariableNames','Delimiter',',')

%%
%get the test matrix for all test, not only math test
data_all_test = [];
VariableNames = [];
Variable_class = []; 
curr_VariableNames = 1;
curr_class = 1;
for curr_task = 1:length(task_data)
    curr_task_fields_list = fields(task_data(curr_task).sesT1);
    for curr_task_field = 1:length(curr_task_fields_list)
        task_field = curr_task_fields_list{curr_task_field};
        if ~contains(task_field,'participant_id')
            if ischar(eval(['task_data(curr_task).sesT1.' task_field '']))
                temp = [];
                temp = eval(['get_number(task_data(curr_task).sesT1.' task_field ')']);
                data_all_test = [data_all_test temp(:)];
            else
                temp = [];
                temp = eval(['task_data(curr_task).sesT1.' task_field ';']);
                data_all_test = [data_all_test temp(:)];
            end
            Variable_class(curr_VariableNames) = curr_class;
            VariableNames{curr_VariableNames} = task_field;
            curr_VariableNames = curr_VariableNames + 1;
        end
    end
    curr_class = curr_class + 1;
end
all_mapID = [1:size(data_all_test,1)];

char_VariableNames = char(VariableNames);
char_VariableNames = char_VariableNames(logical(nansum(data_all_test)),:);
Variable_class = Variable_class(logical(nansum(data_all_test)));
no_nan_sum_all_test = data_all_test(:,logical(nansum(data_all_test)));

% exclude test with too many missing data, standard is boxplot outliers
num_of_nan = sum(isnan(no_nan_sum_all_test));
boxplot(num_of_nan')
remain_ID = num_of_nan<8;
Variable_class = Variable_class(remain_ID);
char_VariableNames = char_VariableNames(remain_ID,:);
char_VariableNames = cellstr(char_VariableNames);

no_nan_sum_all_test = no_nan_sum_all_test(:,remain_ID);

all_mapID = all_mapID(all(~isnan(no_nan_sum_all_test),2));
no_nan_sum_all_test = no_nan_sum_all_test(all(~isnan(no_nan_sum_all_test),2),:);
table = [];
table = array2table(no_nan_sum_all_test,'VariableNames',char_VariableNames);
% writetable(table,'/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan','Delimiter',',')
save('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header','table','all_mapID','Variable_class','char_VariableNames')

table_VariableNames = cell2table(char_VariableNames');
% writetable(table_VariableNames,'/home/kailong/Scheinost-Lab/math/data/VariableNames','Delimiter',',')

%%
clear all;
load('/home/kailong/Scheinost-Lab/math/data/math_test_no_norm_no_nan_header')
% load('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header')
no_norm_no_nan_sum_math_matrix = table2array(math_table);
norm_no_nan_sum_math_matrix = normalize(table2array(math_table),1);

%pick only the test selected here 
% remainTestList = ["CMAT_BasicCalc_Comp_Quotient","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA-2_Attitudes_StS"];TOMA0x2D2_Attitudes_StS
% allRemainTestList = ["ADHD_Total_%","AWMA-S_VerbalSTM_StS","AWMA-S_VerbalWM_StS","AWMA-S_VisuoSpatialSTM_StS","AWMA-S_VisuoSpatialWM_StS","CMAT_BasicCalc_Comp_Quotient","CTOPP_EL_StS","CTOPP_BW_StS","CTOPP_PhonAwareness_Comp","CTOPP_RapidNaming_Comp","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA-2_Attitudes_StS","TOWRE_Total_StS","WASI_Vocab_T-Score","WASI_BD_T-Score","WASI_Sim_T-Score","WASI_MR_T-Score","WASI_VIQ","WASI_PIQ","WASI_FSIQ","WJ-III_WordID_StS","WJ-III_WA_StS","WJ-III_PassComp_StS","WJ-III_MathFluency_StS","WJ-III_SpatialRelations_StS","WJ-III_BRS_StS"];
remainTestList = ["CMAT_BasicCalc_Comp_Quotient","KeyMath_Numeration_ScS","KeyMath_Measurement_ScS","KeyMath_ProblemSolving_ScS","TOMA0x2D2_Attitudes_StS"];
remainTestID = cellfun(@(x) sum(contains(remainTestList,x)),char_math_VariableNames);
remainTestID = find(remainTestID==1);
char_math_VariableNames = char(char_math_VariableNames);
remainTest_math_VariableNames = char_math_VariableNames(remainTestID,:);
remainTest_math_VariableNames = cellstr(remainTest_math_VariableNames);
remainTest_no_norm_no_nan_sum_math_matrix = no_norm_no_nan_sum_math_matrix(:,remainTestID);
remainTest_norm_no_nan_sum_math_matrix = normalize(remainTest_no_norm_no_nan_sum_math_matrix,1);
% remainTest_math_VariableNames = char_math_VariableNames{find(remainTestID==1)};

%check the correlation between different tests
[cr,~] = xcorr(remainTest_norm_no_nan_sum_math_matrix,0,'coeff');
cr = reshape(cr,[sqrt(length(cr)),sqrt(length(cr))]);
cr = cr - diag(diag(cr));
figure
imagesc(cr)

tril_cr = tril(cr);
tril_cr = tril_cr(:);
temp = tril_cr(tril_cr~=0);
figure
boxplot(temp(:))
figure
hist(temp(:),25)

[noHighCor_math_table,noHighCor_math_table_VariableNames] = excludeTooHighCorr(cr,math_Variable_class,char_math_VariableNames,norm_no_nan_sum_math_matrix,no_norm_no_nan_sum_math_matrix);
writetable(noHighCor_math_table,'/home/kailong/Scheinost-Lab/math/data/noHighCor_math_no_norm_no_nan_header','Delimiter',',')
writetable(noHighCor_math_table_VariableNames,'/home/kailong/Scheinost-Lab/math/data/noHighCor_math_VariableNames','Delimiter',',')

%{
clear all;
load('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header')
no_norm_no_nan_sum_all_test = table2array(table);
norm_no_nan_sum_all_test = normalize(table2array(table),1);

%check the correlation between different tests
[cr,~] = xcorr(norm_no_nan_sum_all_test,0,'coeff');
cr = reshape(cr,[sqrt(length(cr)),sqrt(length(cr))]);
cr = cr - diag(diag(cr));

[noHighCor_table,noHighCor_table_VariableNames] = excludeTooHighCorr(cr,Variable_class,char_VariableNames,norm_no_nan_sum_all_test,no_norm_no_nan_sum_all_test);

%}

%pca
% [pca_coefficient,score,latent] = pca(norm_no_nan_sum_math_matrix,'algorithm','als');
[pca_coefficient,score,latent] = pca(norm_no_nan_sum_math_matrix);
%for pruning PC
varianceThreshold = 0.5;
minComponent = 1;
while sum(latent(1:minComponent))/sum(latent(:)) < varianceThreshold %threshold for variance to determine minimum amount of principal components
    minComponent = minComponent+1;
end

summary = rCPM_pca(score,math_mapID);
savedir = '/home/kailong/Scheinost-Lab/math/plot/rCPM/pca_no_nan/pearson/pca/math/';
save([savedir 'performance of pc'],'summary')

figure;
boxplot(cell2mat(summary.q_s));
figure;
boxplot(cell2mat(summary.r_pearson));
figure;
boxplot(cell2mat(summary.r_rank));

pca_for_all_test

% % this is potentially useful for CCA analysis where multiple connectome as
% % well as multiple behavior test can be feed into CPM
% use_scores(:,:) = score(:,1:minComponent);
% pca_coefficient(:,1);

%%

%factor analysis
a = [];
for num_common_factor = 1:10
    [a{num_common_factor}.lambda,a{num_common_factor}.psi,a{num_common_factor}.T,...
        a{num_common_factor}.stats,a{num_common_factor}.F] = ...
        factoran(norm_no_nan_sum_math_matrix,num_common_factor);
%     if stats.p < 0.05
%         break;
%     end
end

% how to derive common_factor? The following should be wrong
common_factor = repmat(lambda',[size(norm_no_nan_sum_math_matrix,1),1]).*norm_no_nan_sum_math_matrix;
common_factor = sum(common_factor,2);
figure;plot([1:17],lambda,'r.');hold on;plot([1:17],1-psi,'b.');plot([1:17],F,'g.');


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

%%
%old pca: only using the one of the connectome 
% clear sum_performance
% savedir = '/home/kailong/Scheinost-Lab/math/plot/pca_no_nan/pearson/pca/';
% if ~isdir(savedir); mkdir(savedir); end
% for pc = 1:size(score,2)
%     all_behav = score(:,pc);
%     temp = [];
%     for i=1:10
%         [y_predict,performance] = cpm_main(all_mats(:,:,math_mapID),all_behav,corr_type,LinearFlag,'pthresh',p,'kfolds',kfold);
%         temp = [temp;performance];
%     end
%     close all;
%     figure;
%     boxplot(temp)
%     sum_performance{pc} = temp; temp = [];
%     savefig([savedir 'performance of pc ' num2str(pc)])
% end
% save([savedir 'performance of pc'],'sum_performance')
% 
% for pc = 1:17
%     sum_r(:,pc) = sum_performance{pc}(:,1);
%     sum_p(:,pc) = sum_performance{pc}(:,2);
% end
% figure;
% boxplot(sum_r)
% figure;
% boxplot(sum_p)
% sig_ID = and(ttest(sum_p-0.05),mean(sum_p)<0.05);
% figure
% boxplot(sum_r(:,sig_ID))
% figure
% boxplot(sum_p(:,sig_ID))

%%
% % [lambda,psi,T,stats,F] = factoran(norm_no_nan_sum_math_matrix,num_common_factor);
% all_behav = common_factor;
% all_behav = F;
% sum_performance=[];
% for i=1:10
%     [y_predict,performance] = cpm_main(all_mats(:,:,math_mapID),all_behav,corr_type,LinearFlag,'pthresh',p,'kfolds',kfold);
%     sum_performance = [sum_performance;performance];
% end
% close all;
% figure;
% boxplot(sum_performance)
% savefig('/home/kailong/Scheinost-Lab/math/plot/pearson/factor analysis/performance of predicted scores')
