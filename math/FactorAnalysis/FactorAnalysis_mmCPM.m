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
AllQs(AllQs==0) = nan;
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
