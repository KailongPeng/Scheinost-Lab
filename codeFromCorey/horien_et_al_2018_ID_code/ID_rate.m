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

x_ses1_original(:,:,motion_ses1,:) = [];
x_ses2_original(:,:,motion_ses2,:) = [];

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

% save(['/home/kailong/Scheinost-Lab/codeFromCorey/horien_et_al_2018_ID_code/ID_code/data/output' ...
%     datestr(datetime('now'))],'sum')
save(['/home/kailong/Scheinost-Lab/codeFromCorey/horien_et_al_2018_ID_code/ID_code/data/output_NoHighMotion2'],'sum')
%%
addpath('/home/kailong/Scheinost-Lab/')
load(['/home/kailong/Scheinost-Lab/codeFromCorey/horien_et_al_2018_ID_code/ID_code/data/output' ...
    '23-Apr-2019 04:40:14.mat'],'sum')


ID_rate1 = cell2mat(kailong_extractfield(sum,'ID_rate1'));
ID_rate2 = cell2mat(kailong_extractfield(sum,'ID_rate2'));

original_ID_rate = [ID_rate1;ID_rate2];

PermutationRate = kailong_cell2mat(kailong_extractfield(sum,'PermutationRate'));
Permute_PermutationRate = permute(PermutationRate,[1 3 2]);
mean_PermutationRate = permute(mean(Permute_PermutationRate),[3 2 1]);

BoostrapRate = kailong_cell2mat(kailong_extractfield(sum,'BoostrapRate'));size(BoostrapRate)
Permute_BoostrapRate = permute(BoostrapRate,[1 3 2]);
mean_BoostrapRate = permute(mean(Permute_BoostrapRate),[3 2 1]);


% mean_PermutationRate = squeeze(mean(PermutationRate));


% PermutationRate_t = reshape(mean(PermutationRate),[2,7]);size(PermutationRate);t = mean(PermutationRate);
% PermutationRate = reshape(PermutationRate,[1000,7,2]);
% BoostrapRate = kailong_cell2mat(kailong_extractfield(sum,'BoostrapRate'));size(BoostrapRate)
% mean_BoostrapRate = squeeze(mean(BoostrapRate));

BoostrapRate = reshape(BoostrapRate,[1000,7,2]);%figure;boxplot(BoostrapRate(:,:,1));
% BoostrapRate = reshape(mean(BoostrapRate),[2,7]);

figure;
subplot(2,3,1);
hold on;
plot(original_ID_rate(1,:),'ro')
plot(original_ID_rate(2,:),'bo')
title('original ID rate')
ylim([0 1])
xticks([1:7])
xticklabels(x_label);xtickangle(45);

subplot(2,3,2);
hold on;
plot(mean_PermutationRate(1,:),'ro')
plot(mean_PermutationRate(2,:),'bo')
title('mean Permutation Rate')
ylim([0 1])
xticks([1:7])
xticklabels(x_label);xtickangle(45);

subplot(2,3,3);hold on;
plot(mean_BoostrapRate(1,:),'ro')
plot(mean_BoostrapRate(2,:),'bo')
title('mean Boostrap Rate')
ylim([0 1])
xticks([1:7])
xticklabels(x_label);xtickangle(45);

subplot(2,3,5);hold on;
boxplot(Permute_PermutationRate(:,:,1),'Colors','b')
boxplot(Permute_PermutationRate(:,:,2),'Colors','g')
title('Permutation Rate')
ylim([0 1])
xticks([1:7])
xticklabels(x_label);xtickangle(45);

subplot(2,3,6);hold on;
boxplot(Permute_BoostrapRate(:,:,1),'Colors','b')
boxplot(Permute_BoostrapRate(:,:,2),'Colors','g')
title('Boostrap Rate')
ylim([0 1])
x_label = ["Mult" "Mult" "Num" "Num" "Rhyming" "Sub" "Sub"];
xticklabels(x_label);xtickangle(45);



saveas(gcf,'summary.png')
% saveas(gcf,'BoostrapRate.png')
% close;

%     x_ses1 = reshape(triu(x_ses1),[],size(x_ses1,3));
%     x_ses2 = reshape(triu(x_ses2),[],size(x_ses2,3));

%
% t=reshape(triu(x_ses1(:,:,1)),[],1);
% t=triu(x_ses1(:,:,1));

