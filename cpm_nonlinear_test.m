clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab/'));

%pathname = '/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST2_LR/matrices/';

kailongLog = [];
curr_log = 1;
task_list = {'REST_LR' 'REST_RL' 'REST2_LR' 'REST2_RL'};
method_list = {'full_no_max' 'partial_correlation' 'full_correlation' };
id = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_id.mat');
id = id.all_id;
id = sort(id);
GSR_list = {'GSR' 'NOGSR'};
save_path = '/home/kailong/Desktop/results_matrix_268_110817/';
if ~isdir(save_path)
    mkdir(save_path);
end
for curr_method = 1:length(method_list)
    method = method_list{curr_method};
    
    missingnodes_hcp = load('/home/kailong/Scheinost-Lab/missingNodes.mat');
    missingnodes_hcp = missingnodes_hcp.missingnodes_hcp;
    behavior = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_behav.mat');
    behavior = behavior.all_behav;
    for curr_GSR = 1:size(GSR_list,2)
        GSR = GSR_list{curr_GSR};
        for curr_task = 1:size(task_list,2)
            clear task task_matrice sum_shift
            task = task_list{curr_task};
            task_matrice = ['MxMxN_matrix_' task];
            savefolder = []; 
            switch method
                case 'partial_correlation'
                    savefolder = [save_path GSR '/partial_correlation/'];
                case 'full_correlation'
                    savefolder = [save_path GSR '/'];
                case 'full_no_max'
                    savefolder = [save_path GSR '/full_no_max/'];
                case 'icov001'
                    savefolder = [save_path GSR '/icov001/'];
                case 'icov005'
                    savefolder = [save_path GSR '/icov005/'];
                otherwise
                    error('method wrong')
            end
            
            output_path = [save_path 'cpm_nonlinear_p_05_k_max/' GSR '/' method '/'];
            output_file = [output_path task_matrice '_cpm_nonlinear'];
            if ~isdir(output_path)
                mkdir(output_path);
            end

            lock_file = [output_file '_lock.mat'];        
            if exist([output_file '.mat'])
                fprintf([task_matrice ' exist\n'])
                continue
            end
            if exist(lock_file)
                fprintf('occupied\n');
                continue;
            else
                save(lock_file,'lock_file');
            end

            clear MxMxN_matrix_temp sum_shift index all_index fileList y summary_r summary_p
            load([savefolder task_matrice]);

            eval(['MxMxN_matrix_temp = ' task_matrice ';']);
            switch method
                case 'partial_correlation'
                    MxMxN_matrix_temp = atanh(MxMxN_matrix_temp(:,:,~isnan(MxMxN_matrix_temp(1,1,:))));
                case 'full_correlation'
                    MxMxN_matrix_temp = MxMxN_matrix_temp(:,:,~isnan(MxMxN_matrix_temp(1,1,:)));
                case 'full_no_max'
                    MxMxN_matrix_temp = MxMxN_matrix_temp(:,:,~isnan(MxMxN_matrix_temp(1,1,:)));
                case 'icov001'
                    MxMxN_matrix_temp = MxMxN_matrix_temp(:,:,~isnan(MxMxN_matrix_temp(1,1,:)));
                case 'icov005'
                    MxMxN_matrix_temp = MxMxN_matrix_temp(:,:,~isnan(MxMxN_matrix_temp(1,1,:)));
                otherwise
                    error('method wrong')
            end
            
            curr_loop = 1;
%             summary_r = nan(length(p_value_list),length(kfold_list));
%             summary_p = nan(length(p_value_list),length(kfold_list));
            sum_r = [];
            sum_p = [];
            tic
%             for curr_p = 1:size(p_value_list,2)
%                 p = [];
%                 p = p_value_list(curr_p);
%                 for curr_kfold = 1:size(kfold_list,2)
%                     kfold = [];
%                     kfold = kfold_list(curr_kfold);
%                     y{curr_loop}.p = p;
%                     y{curr_loop}.kfold = kfold;
%                     corr_type = 'spearman';
%                     LinearFlag = 1;
%                     try
%                         [y{curr_loop}.y_predict, y{curr_loop}.performance] = cpm_main(MxMxN_matrix_temp,behavior(map_ID),corr_type,LinearFlag,'pthresh',p,'kfolds',kfold);
%                         summary_r(curr_p,curr_kfold) = y{curr_loop}.performance(1);
%                         summary_p(curr_p,curr_kfold) = y{curr_loop}.performance(2);    
%                     end
%                     curr_loop = curr_loop + 1;
%                 end
%             end
            corr_type_List = {'Pearson' 'spearman'};
            LinearFlagList = [1,0];
            curr_loop = 1;
            p = 0.05;
            kfold = size(MxMxN_matrix_temp,3);
%             kfold = 2;
            for curr_corr_type = 1:length(corr_type_List)
                corr_type = corr_type_List{curr_corr_type};
                for curr_LinearFlag = 1:length(LinearFlagList)
                    LinearFlag = LinearFlagList(curr_LinearFlag);
                    [y{curr_loop}.y_predict, y{curr_loop}.performance] = cpm_main(MxMxN_matrix_temp,behavior(map_ID),corr_type,LinearFlag,'pthresh',p,'kfolds',kfold);
                    sum_r(curr_corr_type,curr_LinearFlag) = y{curr_loop}.performance(1);
                    sum_p(curr_corr_type,curr_LinearFlag) = y{curr_loop}.performance(2);
                    curr_loop = curr_loop + 1;
                end
            end
            try
                save(output_file,'y','sum_r','sum_p')
            catch
                warning('unsuccessful fit!')
                kailongLog{curr_log} = [savefolder task_matrice '-unsuccessful fit!' ];
                curr_log = curr_log + 1;
            end
            delete(lock_file)
            toc
        end
    end

%     for curr_log = 1:size(kailongLog,2)
%         fprintf(num2str(kailongLog{curr_log}));
%     end
%     kailongLog = [];
%     curr_log = 1;
end

