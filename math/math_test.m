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

% task_math(1).sesT1.participant_id_num = str2num(task_math(1).sesT1.participant_id(end-3,end));
% temp = task_math(3).sesT1.TOMA0x2D2_Attitudes_StS;
% temp2 = get_number(temp);
sum_math_matrix = [];
math_VariableNames = [];
curr_VariableNames = 1;
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
            math_VariableNames{curr_VariableNames} = task_field;
            curr_VariableNames = curr_VariableNames + 1;
        end
    end
end
char_math_VariableNames = char(math_VariableNames);
char_math_VariableNames = char_math_VariableNames(logical(nansum(sum_math_matrix)),:);
char_math_VariableNames = cellstr(char_math_VariableNames);
no_nan_sum_math_matrix = sum_math_matrix(:,logical(nansum(sum_math_matrix)));

no_nan_sum_math_matrix = no_nan_sum_math_matrix(all(~isnan(no_nan_sum_math_matrix),2),:);
norm_no_nan_sum_math_matrix = normalize(no_nan_sum_math_matrix,1);
math_table = [];
math_table = array2table(norm_no_nan_sum_math_matrix,'VariableNames',char_math_VariableNames);
% writetable(math_table,'/home/kailong/Scheinost-Lab/math/data/math_test_norm_no_nan_header','Delimiter',',')

no_norm_no_nan_sum_math_matrix = no_nan_sum_math_matrix;
math_table = [];
math_table = array2table(no_norm_no_nan_sum_math_matrix,'VariableNames',char_math_VariableNames);
% writetable(math_table,'/home/kailong/Scheinost-Lab/math/data/math_test_no_norm_no_nan_header','Delimiter',',')

norm_sum_math_matrix = sum_math_matrix(:,logical(nansum(sum_math_matrix)));
norm_sum_math_matrix = normalize(norm_sum_math_matrix,1);
norm_math_table = [];
norm_math_table = array2table(norm_sum_math_matrix,'VariableNames',char_math_VariableNames);
% writetable(norm_math_table,'/home/kailong/Scheinost-Lab/math/data/all_test_norm','Delimiter',',')

table_math_VariableNames = cell2table(char_math_VariableNames');
% writetable(table_math_VariableNames,'/home/kailong/Scheinost-Lab/math/data/math_VariableNames','Delimiter',',')


% save('/home/kailong/Scheinost-Lab/math/data/math_test','sum_matrix')
%%
%get the test matrix for all test, not only math test
data_all_test = [];
VariableNames = [];
curr_VariableNames = 1;
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
            VariableNames{curr_VariableNames} = task_field;
            curr_VariableNames = curr_VariableNames + 1;
        end
    end
end

normalization_flag = 0;clear 
no_nan_flag = 1;
if normalization_flag == 0
    if no_nan_flag == 0
        char_VariableNames = char(VariableNames);
        char_VariableNames = char_VariableNames(logical(nansum(data_all_test)),:);
        char_VariableNames = cellstr(char_VariableNames);
        no_nan_sum_all_test = data_all_test(:,logical(nansum(data_all_test)));
        table = [];
        table = array2table(no_nan_sum_all_test,'VariableNames',char_VariableNames);
        writetable(table,'/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_nan','Delimiter',',')
    else
        % exclude test with too many missing data, standard is boxplot outliers
        num_of_nan = sum(isnan(no_nan_sum_all_test));
        boxplot(num_of_nan')
        remain_ID = num_of_nan<8;
        char_VariableNames = char(VariableNames);
        char_VariableNames = char_VariableNames(remain_ID);
        char_VariableNames = cellstr(char_VariableNames);
        no_nan_sum_all_test = no_nan_sum_all_test(:,remain_ID);
        no_nan_sum_all_test = no_nan_sum_all_test(all(~isnan(no_nan_sum_all_test),2),:);
        table = [];
        table = array2table(no_nan_sum_all_test,'VariableNames',char_VariableNames);
%         writetable(table,'/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan','Delimiter',',')
        table_VariableNames = cell2table(char_VariableNames');
%         writetable(table_VariableNames,'/home/kailong/Scheinost-Lab/math/data/VariableNames','Delimiter',',')
    end
else
    no_nan_sum_all_test = no_nan_sum_all_test(all(~isnan(no_nan_sum_all_test),2),:);
    % sum(sum(isnan(no_nan_sum_all_test)))
    norm_no_nan_sum_all_test = normalize(no_nan_sum_all_test,1);
    table = [];
    table = array2table(norm_no_nan_sum_all_test,'VariableNames',char_VariableNames);
    save('/home/kailong/Scheinost-Lab/math/data/all_test_norm_no_nan','table')
    writetable(table,'/home/kailong/Scheinost-Lab/math/data/all_test_norm_no_nan','Delimiter',',')
    
    norm_sum_all_test = data_all_test(:,logical(nansum(data_all_test)));
    norm_sum_all_test = normalize(norm_sum_all_test,1);
    norm_table = [];
    norm_table = array2table(norm_sum_all_test,'VariableNames',char_VariableNames);
    writetable(norm_table,'/home/kailong/Scheinost-Lab/math/data/all_test_norm','Delimiter',',')
    
end
% csvwrite('/home/kailong/Scheinost-Lab/math/data/all_test_norm_no_nan',table)

%%
load('/home/kailong/Scheinost-Lab/math/data/math_test')
sum_matrix_norm = normalize(sum_math_matrix,1);
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

%%
clear sum_performance
savedir = '/home/kailong/Scheinost-Lab/math/plot/pca_no_nan/pearson/pca/';
if ~isdir(savedir); mkdir(savedir); end;
for pc = 1:size(score,2)
    all_behav = score(:,pc);
    temp = [];
    for i=1:10
        [y_predict,performance] = cpm_main(all_mats(:,:,mapID),all_behav,corr_type,LinearFlag,'pthresh',p,'kfolds',kfold);
        temp = [temp;performance];
    end
    close all;
    figure;
    boxplot(temp)
    sum_performance{pc} = temp; temp = [];
    savefig([savedir 'performance of pc ' num2str(pc)])
end
save([savedir 'performance of pc'],'sum_performance')

for pc = 1:17
    sum_r(:,pc) = sum_performance{pc}(:,1);
    sum_p(:,pc) = sum_performance{pc}(:,2);
end
figure;
boxplot(sum_r)
figure;
boxplot(sum_p)
sig_ID = and(ttest(sum_p-0.05),mean(sum_p)<0.05);
figure
boxplot(sum_r(:,sig_ID))
figure
boxplot(sum_p(:,sig_ID))

%%
% [lambda,psi,T,stats,F] = factoran(sum_matrix_norm_no_nan,num_common_factor);
all_behav = common_factor;
all_behav = F;
sum_performance=[];
for i=1:10
    [y_predict,performance] = cpm_main(all_mats(:,:,mapID),all_behav,corr_type,LinearFlag,'pthresh',p,'kfolds',kfold);
    sum_performance = [sum_performance;performance];
end
close all;
figure;
boxplot(sum_performance)
savefig('/home/kailong/Scheinost-Lab/math/plot/pearson/factor analysis/performance of predicted scores')
