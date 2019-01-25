%Matrices and ROI mean text files (both GSR and noGSR) for each subject 
%can be found in /data1/software/HCP_data/HCP_900_DATA/TASK_[RL/LR]/matrices/###*.txt 
%with ### = subject ID. 


%IDs for the 515 subjects we used in the NComms paper are in 
%/mnt/store4/mri_group/siyuan_data/HCP515/all_id.mat and pmat scores 
%for those subjects can be found in /mnt/store4/mri_group/siyuan_data/HCP515/all_behav.mat. 
%Indices of nodes that are missing in any subjects can be found in 
%/mnt/newchell/47421/NeurodevelopmentalGenomics/abby/CPMPaper/hcp515_noBadNodes/missingNodes.mat. 

clear all;close all;clc;

%pathname = '/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST2_LR/matrices/';
pathname = '/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/EMOTION_LR/matrices/';


%cd(pathname);
a = dir([pathname '*_GSR_*']);
filename = [pathname '100206_REST2_LR_GSR_matrix.txt'];
connectivity_matrices = readtable(pathname);
connectivity_matrices(:,end) = [];
connectivity_matrices = table2array(connectivity_matrices);
id = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_id.mat');
id = id.all_id;
behavior = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_behav.mat');
behavior = behavior.all_behav;
% notworking = load('/mnt/newchell/47421/NeurodevelopmentalGenomics/abby/CPMPaper/hcp515_noBadNodes/missingNodes.mat');


[y_predict, performance] = cpm_main(connectivity_matrices,behavior,'pthresh',0.05,'kfolds',2)

