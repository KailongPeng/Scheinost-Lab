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
fileList = kailong_extractfield(a,'name');
id = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_id.mat');
id = id.all_id;
for curr_id = 1:size(id,1) % find the filenameID that is among the id list
    Index{curr_id} = find(~cellfun(@isempty,strfind(fileList, num2str(id(curr_id)))));
end
for curr_id = 1:size(id,1) % find the filenameID that is among the id list
    Index{curr_id} = find(contains(fileList,num2str(id(curr_id)))==1);
end



find(strcmp(

filename = [pathname '100206_REST2_LR_GSR_matrix.txt'];
connectivity_matrices = readtable(pathname);
connectivity_matrices(:,end) = [];
connectivity_matrices = table2array(connectivity_matrices);
behavior = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_behav.mat');
behavior = behavior.all_behav;
% notworking = load('/mnt/newchell/47421/NeurodevelopmentalGenomics/abby/CPMPaper/hcp515_noBadNodes/missingNodes.mat');


%{
addpath('/home/kailong/Scheinost-Lab')
addpath('/home/kailong/CPM/matlab')
load ('/home/kailong/Desktop/mega_sample_200_female_subjs.mat')


[y_predict, performance] = cpm_main(MxMxN_matrix_FEMALE_HC_200,all_behav_FEMALE_200,'pthresh',0.05,'kfolds',2);
1;
%}
