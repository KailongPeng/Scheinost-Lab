%calculate ID_rate
clear;close all;clc;
restoredefaultpath;
addpath('/home/kailong/Scheinost-Lab/codeFromCorey/horien_et_al_2018_ID_code/ID_code')
x_ses1_original = load('/home/kailong/Scheinost-Lab/math/data/TwoSession/all_mats_ses1.mat');x_ses1_original = x_ses1_original.all_mats;%268node*268node*132sub*7task
x_ses2_original = load('/home/kailong/Scheinost-Lab/math/data/TwoSession/all_mats_ses2.mat');x_ses2_original = x_ses2_original.all_mats;%268node*268node*132sub*7task

incomp_id_ses1 = load('/home/kailong/Scheinost-Lab/math/data/TwoSession/incomp_id_ses1.mat');incomp_id_ses1 = incomp_id_ses1.incomp_id;
incomp_id_ses2 = load('/home/kailong/Scheinost-Lab/math/data/TwoSession/incomp_id_ses2.mat');incomp_id_ses2 = incomp_id_ses2.incomp_id;
motion_ses1 = load('/home/kailong/Scheinost-Lab/math/data/TwoSession/motion_ses1.mat');motion_ses1 = motion_ses1.bad_id;
motion_ses2 = load('/home/kailong/Scheinost-Lab/math/data/TwoSession/motion_ses2.mat');motion_ses2 = motion_ses2.bad_id;
sum = [];
for curr_task = 1:size(x_ses1_original,4)
    fprintf('%d\n',curr_task);
    tic;
    x_ses1 = []; x_ses2 = [];t_ses1 = [];t_ses2 = [];
    x_ses1 = x_ses1_original(:,:,:,curr_task);
    x_ses2 = x_ses2_original(:,:,:,curr_task);
    for curr_sub = 1:size(x_ses1,3)
        x_ses1(:,:,curr_sub) = x_ses1(:,:,curr_sub) - diag(diag(x_ses1(:,:,curr_sub)));
        t_ses1{curr_sub} = reshape(triu(x_ses1(:,:,curr_sub)),[],1);
        %     t_ses1{curr_sub}(t_ses1{curr_sub}==0) = [];
        
        
        x_ses2(:,:,curr_sub) = x_ses2(:,:,curr_sub) - diag(diag(x_ses2(:,:,curr_sub)));
        t_ses2{curr_sub} = reshape(triu(x_ses2(:,:,curr_sub)),[],1);
        %     t_ses2{curr_sub}(t_ses2{curr_sub}==0) = [];
    end
    x_ses1 = cell2mat(t_ses1);
    x_ses2 = cell2mat(t_ses2);

    all_se1 = x_ses1;
    all_se2 = x_ses2;
    
    %ID_rate
    t = [];
    t = ID_example(all_se1,all_se2);  
    sum{curr_task}.ID_rate1 = t.rate1;
    sum{curr_task}.ID_rate2 = t.rate2;
    %permutation_test
    sum{curr_task}.PermutationRate = permutation_test_example(x_ses1,x_ses2);
    %bootstrap_example
    sum{curr_task}.BoostrapRate = bootstrap_example(x_ses1,x_ses2);
    toc
end

save(['/home/kailong/Scheinost-Lab/codeFromCorey/horien_et_al_2018_ID_code/ID_code/data/output' ...
    datestr(datetime('now'))],'sum')
%     x_ses1 = reshape(triu(x_ses1),[],size(x_ses1,3));
%     x_ses2 = reshape(triu(x_ses2),[],size(x_ses2,3));

%
% t=reshape(triu(x_ses1(:,:,1)),[],1);
% t=triu(x_ses1(:,:,1));

