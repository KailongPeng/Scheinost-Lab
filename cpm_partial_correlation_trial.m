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

kailongLog = [];
curr_log = 1;
% task_list = {'REST_LR' 'REST_RL' 'REST2_LR' 'REST2_RL' ...
%     'EMOTION_LR' 'EMOTION_RL' ...
%     'MOTOR_LR' 'MOTOR_RL'...
%     'SOCIAL_LR' 'SOCIAL_RL'...
%     'WM_LR' 'WM_RL'...
%     'GAMBLING_LR' 'GAMBLING_RL'...
%     'LANGUAGE_LR' 'LANGUAGE_RL'...
%     'RELATIONAL_LR' 'RELATIONAL_RL'...
%     }
task_list = {'REST_LR' 'REST_RL' 'REST2_LR' 'REST2_RL'};
method_list = {'partial_correlation' 'full_correlation' 'full_no_max' 'icov001' 'icov005' };
id = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_id.mat');
id = id.all_id;
id = sort(id);
GSR_list = {'GSR' 'NOGSR'};
for curr_method = 1:length(method_list)
    method = method_list{curr_method};
    for curr_GSR = 1:size(GSR_list,2)
        GSR = GSR_list{curr_GSR};
        for curr_task = 1:size(task_list,2)
            clear task task_matrice
            task = task_list{curr_task};
            task_matrice = ['MxMxN_matrix_' task];
            savefolder = []; 
            switch method
                case 'partial_correlation'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/partial_correlation/'];
                case 'full_correlation'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/'];
                case 'full_no_max'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/full_no_max/'];
                case 'icov001'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/icov001/'];
                case 'icov005'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/icov005/'];
                otherwise
                    error('method wrong')
            end
            
            if ~isdir(savefolder)
                mkdir(savefolder)
            end
            if exist([savefolder task_matrice '.mat'])
%                 load([savefolder task_matrice '.mat'])
%                 clear temp
%                 eval(['temp = ' task_matrice ';']);
%                 temp = temp(~isnan(temp));
%                 if ~isempty(temp)
                    fprintf('exist\n')
                    continue
%                 end
%                 clear temp
            end
            lock_file = [savefolder task_matrice '_lock.mat'];        
            if exist(lock_file)
                fprintf('occupied\n');
                continue
            else
                save(lock_file,'lock_file');
            end
            pathname = ['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/' task '/matrices/'];
            cd(pathname);
            clear a fileList;
            a = dir([pathname '*_' GSR '_roimean*']);
            fileList = kailong_extractfield(a,'name');
            clear a;

            clear index all_index map_ID
            % for curr_id = 1:size(id,1) % find the filenameID that is among the id list
            %     Index{curr_id} = find(~cellfun(@isempty,strfind(fileList, num2str(id(curr_id)))));
            % end
            for curr_id = 1:size(id,1) % find the filenameID that is among the id list
                if sum(contains(fileList,num2str(id(curr_id)))) ~= 0
                    all_index{curr_id} = find(contains(fileList,num2str(id(curr_id)))==1);
                else
                    all_index{curr_id} = nan;
                    warning('missing subjects')
                    kailongLog{curr_log} = ['missing subjects\n'];
                    curr_log = curr_log + 1; 
                end
            end
            all_index = cell2mat(all_index);
            map_ID = [1:size(all_index,2)];

            map_ID = map_ID(~isnan(all_index));
            index = all_index(~isnan(all_index));
            new_file_list = [];
            for curr_sub = 1:size(index,2)
                new_file_list{curr_sub} = fileList{index(curr_sub)};
            end

            clear  filename shift cr_max MxMxN_matrix_temp sum_shift  %curr_data cr lgs corrected_cr_max
            Log = [];
            eval(['clear ' task_matrice])

    %         if exist([savefolder task_matrice '.mat'])
    %             fprintf('exist\n')
    %             save([savefolder task_matrice],'index','all_index','fileList','map_ID','-append');
    %             continue
    %         end
            switch method
                case 'partial_correlation'

                    parfor (curr_sub = 1:size(new_file_list,2),4)
                        curr_sub
                        filename = [pathname new_file_list{curr_sub}];
                        try % discard subjects with missing nodes
                            [rho,atanh_rho] = create_partial_correlation_matrix_function(filename);
                            MxMxN_matrix_temp(:,:,curr_sub) = rho;
                            MxMxN_matrix_temp_norm(:,:,curr_sub) = atanh_rho;
                            filename = [];
                            rho = []; 
                            atanh_rho = [];
                        catch
                            map_ID(curr_sub) = nan;
                            MxMxN_matrix_temp(:,:,curr_sub) = nan;
                            MxMxN_matrix_temp_norm(:,:,curr_sub) = nan;

                            filename = [];
                            rho = []; 
                            atanh_rho = [];
                            warning('jumping subject')
                            Log = [Log 'jumping subject\n'];
                        end
                    end
                    map_ID = map_ID(~isnan(map_ID));
                    kailongLog{curr_log} = Log;
                    curr_log = curr_log + 1; 

                    %check whether empty
                    clear temp
                    temp = MxMxN_matrix_temp;
                    temp = temp(~isnan(temp));
                    if isempty(temp)
                        error('empty matrix, check Log')
                    end
                    clear temp

                    eval([task_matrice ' = MxMxN_matrix_temp;'])
                    eval([task_matrice '_norm = MxMxN_matrix_temp_norm;'])
                    save([savefolder task_matrice],task_matrice,[task_matrice '_norm']);
                    save([savefolder task_matrice],'index','all_index','fileList','map_ID','-append');
                    delete(lock_file)
                case 'full_correlation'

                    parfor curr_sub = 1:size(new_file_list,2)
                        curr_sub
                        filename = [pathname new_file_list{curr_sub}];
                        try % discard subjects with missing nodes
                            [cr_max,shift] = create_max_correlation_matrix(filename);
                            MxMxN_matrix_temp(:,:,curr_sub) = cr_max;
                            filename = [];
                            sum_shift(:,curr_sub) = shift;
                            shift = [];
                            cr_max = []; %curr_data cr lgs corrected_cr_max
                        catch
                            map_ID(curr_sub) = nan;
                            MxMxN_matrix_temp(:,:,curr_sub) = nan;
                            filename = [];
                            sum_shift(:,curr_sub) = nan;
                            shift = [];
                            cr_max = []; %curr_data cr lgs corrected_cr_max
                            warning('jumping subject')
                            Log = [Log 'jumping subject\n'];
                        end
                    end
                    map_ID = map_ID(~isnan(map_ID));
                    kailongLog{curr_log} = Log;
                    curr_log = curr_log + 1;          
                    eval([task_matrice ' = MxMxN_matrix_temp;'])

                    save([savefolder task_matrice],task_matrice);
                    save([savefolder task_matrice],'sum_shift','index','all_index','fileList','map_ID','-append');
                    delete(lock_file)
                case 'full_no_max'
%                     parfor curr_sub = 1:size(new_file_list,2)
                    parfor curr_sub = 1:size(new_file_list,2)
                        curr_sub
                        filename = [pathname new_file_list{curr_sub}];
                        try % discard subjects with missing nodes
                            [atanh_cr,cr] = create_no_max_correlation_matrix(filename);
                            MxMxN_matrix_temp(:,:,curr_sub) = atanh_cr;
                            filename = [];
                        catch
                            map_ID(curr_sub) = nan;
                            MxMxN_matrix_temp(:,:,curr_sub) = nan;
                            filename = [];
                            warning('jumping subject')
                            Log = [Log 'jumping subject\n'];
                        end
                    end
                    map_ID = map_ID(~isnan(map_ID));
                    kailongLog{curr_log} = Log;
                    curr_log = curr_log + 1;          
                    eval([task_matrice ' = MxMxN_matrix_temp;'])

                    save([savefolder task_matrice],task_matrice);
                    save([savefolder task_matrice],'index','all_index','fileList','map_ID','-append');
                    delete(lock_file)  
                case 'icov001'
%                     parfor curr_sub = 1:size(new_file_list,2)
                    parfor curr_sub = 1:size(new_file_list,2)
%                         curr_sub
                        tic
                        filename = [pathname new_file_list{curr_sub}];
                        try % discard subjects with missing nodes
                            [atanh_cr,cr] = create_no_max_correlation_matrix(filename);
                            temp = [];
                            temp = L1precisionBCD(cr,.001);
                            temp = real(temp);
                            temp = temp - diag(diag(temp));
                            temp = atanh(temp);
                            MxMxN_matrix_temp(:,:,curr_sub) = temp;
                            temp = []; 
                            filename = [];
                        catch
                            map_ID(curr_sub) = nan;
                            MxMxN_matrix_temp(:,:,curr_sub) = nan;
                            filename = [];
                            warning('jumping subject')
                            Log = [Log 'jumping subject\n'];
                        end
                        toc
                    end
                    map_ID = map_ID(~isnan(map_ID));
                    kailongLog{curr_log} = Log;
                    curr_log = curr_log + 1;          
                    eval([task_matrice ' = MxMxN_matrix_temp;'])

                    save([savefolder task_matrice],task_matrice);
                    save([savefolder task_matrice],'index','all_index','fileList','map_ID','-append');
                    delete(lock_file)    
                case 'icov005'
                    parfor curr_sub = 1:size(new_file_list,2)
                        curr_sub
                        filename = [pathname new_file_list{curr_sub}];
                        try % discard subjects with missing nodes
                            [atanh_cr,cr] = create_no_max_correlation_matrix(filename);
                            temp = [];
                            temp = L1precisionBCD(cr,.005);
                            temp = real(temp);
                            temp = temp - diag(diag(temp));
                            temp = atanh(temp);
                            MxMxN_matrix_temp(:,:,curr_sub) = temp;
                            temp = []; 
                            filename = [];
                        catch
                            map_ID(curr_sub) = nan;
                            MxMxN_matrix_temp(:,:,curr_sub) = nan;
                            filename = [];
                            warning('jumping subject')
                            Log = [Log 'jumping subject\n'];
                        end
                    end
                    map_ID = map_ID(~isnan(map_ID));
                    kailongLog{curr_log} = Log;
                    curr_log = curr_log + 1;          
                    eval([task_matrice ' = MxMxN_matrix_temp;'])

                    save([savefolder task_matrice],task_matrice);
                    save([savefolder task_matrice],'index','all_index','fileList','map_ID','-append');
                    delete(lock_file)   
                otherwise
                        
            end
        end 
    end
 
    for curr_log = 1:size(kailongLog,2)
        fprintf(num2str(kailongLog{curr_log}));
    end
    format long g
    save(['/home/kailong/Scheinost-Lab/Log/' datestr(datetime) '_Log'],'kailongLog');
    kailongLog = [];
    curr_log = 1;
    % fprintf('%s',num2str(datenum(datetime)))
end
%%

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
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/partial_correlation/'];
                case 'full_correlation'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/'];
                case 'full_no_max'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/full_no_max/'];
                case 'icov001'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/icov001/'];
                case 'icov005'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/icov005/'];
                otherwise
                    error('method wrong')
            end
            
            output_file = [savefolder task_matrice '_atanh_p5'];

            lock_file = [output_file '_lock.mat'];        
%             if exist([output_file '.mat'])
%                 fprintf([task_matrice ' exist\n'])
%                 continue
%             end
%             if exist(lock_file)
%                 fprintf('occupied\n');
%                 continue;
%             else
%                 save(lock_file,'lock_file');
%             end

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
            
            p_value_list = [0.5 0.3 0.05 0.01 0.005 0.001];
            kfold_list = [2,10,size(MxMxN_matrix_temp,3)];
            curr_loop = 1;
            summary_r = nan(length(p_value_list),length(kfold_list));
            summary_p = nan(length(p_value_list),length(kfold_list));

            tic
            for curr_p = 1:size(p_value_list,2)
                p = [];
                p = p_value_list(curr_p);
                for curr_kfold = 1:size(kfold_list,2)
                    kfold = [];
                    kfold = kfold_list(curr_kfold);
                    y{curr_loop}.p = p;
                    y{curr_loop}.kfold = kfold;
                    corr_type = 'spearman';
                    LinearFlag = 1;
                    try
                        [y{curr_loop}.y_predict, y{curr_loop}.performance] = cpm_main(MxMxN_matrix_temp,behavior(map_ID),corr_type,LinearFlag,'pthresh',p,'kfolds',kfold);
                        summary_r(curr_p,curr_kfold) = y{curr_loop}.performance(1);
                        summary_p(curr_p,curr_kfold) = y{curr_loop}.performance(2);    
                    end
                    curr_loop = curr_loop + 1;
                end
            end
            try
                save(output_file,'y','summary_r','summary_p')
            catch
                warning('unsuccessful fit!')
                kailongLog{curr_log} = [savefolder task_matrice '-unsuccessful fit!' ];
                curr_log = curr_log + 1;
            end
            delete(lock_file)
            toc
        end
    end

    for curr_log = 1:size(kailongLog,2)
        fprintf(num2str(kailongLog{curr_log}));
    end
    kailongLog = [];
    curr_log = 1;
end

% notworking = load('/mnt/newchell/47421/NeurodevelopmentalGenomics/abby/CPMPaper/hcp515_noBadNodes/missingNodes.mat');
%{
load (['/home/kailong/Desktop/results_matrix_268_110817/MxMxN_matrix_REST_LR.mat'])


load ('/home/kailong/Desktop/mega_sample_200_female_subjs.mat')
[y_predict, performance] = cpm_main(MxMxN_matrix_FEMALE_HC_200,all_behav_FEMALE_200,'pthresh',0.05,'kfolds',2);
performance
%}


