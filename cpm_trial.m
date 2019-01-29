%Matrices and ROI mean text files (both GSR and noGSR) for each subject 
%can be found in /data1/software/HCP_data/HCP_900_DATA/TASK_[RL/LR]/matrices/###*.txt 
%automatically redirected to '/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA'
%with ### = subject ID. 

%IDs for the 515 subjects we used in the NComms paper are in 
%/mnt/store4/mri_group/siyuan_data/HCP515/all_id.mat and pmat scores
%for those subjects can be found in /mnt/store4/mri_group/siyuan_data/HCP515/all_behav.mat. 

%Indices of nodes that are missing in any subjects can be found in 
%/mnt/newchell/47421/NeurodevelopmentalGenomics/abby/CPMPaper/hcp515_noBadNodes/missingNodes.mat. 

clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab/'));
%pathname = '/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST2_LR/matrices/';

task_list = {'REST_LR' 'REST_RL' 'REST2_LR' 'REST2_RL' ...
    'EMOTION_LR' 'EMOTION_RL' ...
    'MOTOR_LR' 'MOTOR_RL'...
    'SOCIAL_LR' 'SOCIAL_RL'...
    'WM_LR' 'WM_RL'...
    'GAMBLING_LR' 'GAMBLING_RL'...
    'LANGUAGE_LR' 'LANGUAGE_RL'...
    'RELATIONAL_LR' 'RELATIONAL_RL'...
    };
% task_list = {'REST_LR' 'REST_RL' 'REST2_LR' 'REST2_RL'};

id = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_id.mat');
id = id.all_id;
id = sort(id);
% GSR_list = {'GSR' 'NOGSR'};
GSR_list = {'GSR'};

for curr_GSR = 1:size(GSR_list,2)
    GSR = GSR_list{curr_GSR};
    for curr_task = 1:size(task_list,2)
        clear task task_matrice
        task = task_list{curr_task};
        task_matrice = ['MxMxN_matrix_' task];
        savefolder = []; 
        savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/'];
        if exist([savefolder task_matrice '.mat'])
            fprintf('exist\n')
            continue
        end
        pathname = ['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/' task '/matrices/'];
        cd(pathname);
        clear a;
        a = dir([pathname '*_' GSR '_roimean*']);
        fileList = kailong_extractfield(a,'name');
        clear a;

        % for curr_id = 1:size(id,1) % find the filenameID that is among the id list
        %     Index{curr_id} = find(~cellfun(@isempty,strfind(fileList, num2str(id(curr_id)))));
        % end
        for curr_id = 1:size(id,1) % find the filenameID that is among the id list
            index{curr_id} = find(contains(fileList,num2str(id(curr_id)))==1);
        end
        new_file_list = [];
        for curr_sub = 1:size(id,1)
            new_file_list{curr_sub} = fileList{index{curr_sub}};
        end

        clear  filename shift cr_max MxMxN_matrix_temp %curr_data cr lgs corrected_cr_max
        eval(['clear ' task_matrice])
        parfor curr_sub = 1:size(new_file_list,2)
            curr_sub
            filename = [pathname new_file_list{curr_sub}];
            [cr_max,shift] = create_max_correlation_matrix(filename);
            MxMxN_matrix_temp(:,:,curr_sub) = cr_max;
            filename = [];
            shift = [];
            cr_max = []; %curr_data cr lgs corrected_cr_max
        end
        eval([task_matrice ' = MxMxN_matrix_temp;'])
        if ~isdir(savefolder)
            mkdir(savefolder)
        end
        save([savefolder task_matrice],task_matrice);
    end
end
missingnodes_hcp = load('/home/kailong/Scheinost-Lab/missingNodes.mat');
missingnodes_hcp = missingnodes_hcp.missingnodes_hcp;
behavior = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_behav.mat');
behavior = behavior.all_behav;
for curr_GSR = 1:size(GSR_list,2)
    GSR = GSR_list{curr_GSR};
    for curr_task = 1:size(task_list,2)
        clear task task_matrice
        task = task_list{curr_task};
        task_matrice = ['MxMxN_matrix_' task];
        savefolder = []; 
        savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/'];
        load([savefolder task_matrice]);
        clear MxMxN_matrix_temp
        eval(['MxMxN_matrix_temp = ' task_matrice ';']);
        [y_predict, performance] = cpm_main(MxMxN_matrix_temp,behavior,'pthresh',0.05,'kfolds',2);
        save([savefolder task_matrice],'y_predict','performance','-append')
    end
end

% notworking = load('/mnt/newchell/47421/NeurodevelopmentalGenomics/abby/CPMPaper/hcp515_noBadNodes/missingNodes.mat');
%{
load (['/home/kailong/Desktop/results_matrix_268_110817/MxMxN_matrix_REST_LR.mat'])


load ('/home/kailong/Desktop/mega_sample_200_female_subjs.mat')
[y_predict, performance] = cpm_main(MxMxN_matrix_FEMALE_HC_200,all_behav_FEMALE_200,'pthresh',0.05,'kfolds',2);
performance
%}
