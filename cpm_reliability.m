
clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab/'));

kailongLog = [];
curr_log = 1;
task_list = {'REST_LR' 'REST_RL' 'REST2_LR' 'REST2_RL'};
method_list = {'full_no_max' 'full_correlation' 'partial_correlation'};
id = load('/mnt/store4/mri_group/siyuan_data/HCP515/all_id.mat');
id = id.all_id;
id = sort(id);
GSR_list = {'GSR' 'NOGSR'};
% partial = input('partial correlation? y or n\n','s');
% if strcmp(partial,'y')
% if contains(filename,'partial')
%     partial = 'partial_correlation';
% else
%     partial = '';
% end
sub_num = input('number of subjects you want to include in the analysis (max 515)\n');
summary = [];
summary_all = [];
for curr_method = 1:length(method_list)
    method = method_list{curr_method};
    for curr_GSR = 1:size(GSR_list,2)
        GSR = GSR_list{curr_GSR};
%         for curr_task = 1:size(task_list,2)
%             clear task task_matrice
%             task = task_list{curr_task};
%             task_matrice = ['MxMxN_matrix_' task];
        savefolder = []; 
        filename = ['cpm_reliability_original_' method '_' GSR 'sub_num_' num2str(sub_num)];
        switch method
            case 'partial_correlation'
                savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/partial_correlation/'];
            case 'full_correlation'
                savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/'];
            case 'full_no_max'
                savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/full_no_max/'];
            otherwise
                error('savefolder = []');
        end

%         if ~isdir(savefolder)
%             mkdir(savefolder)
%         end
%             if exist([savefolder task_matrice '.mat'])
%                     fprintf('exist\n')
%                     continue
%             end
%             lock_file = [savefolder task_matrice '_lock.mat'];        
%             if exist(lock_file)
%                 fprintf('occupied\n');
%                 continue
%             else
%                 save(lock_file,'lock_file');
%             end
        thisFolder = savefolder;
        thisPattern = 'cpm_reliability';
        [data,ftbl] = load_reliability_data_for_CpmReliability(thisFolder, thisPattern,filename,task_list,sub_num);
%             for i = 1:length(task_list)
%                 theseFiles{i} = [thisFolder 'MxMxN_matrix_' task_list{i} '.mat'];
%             end
        correctiontype = 'none';                                                                  
        [icc_summary,var,stats,sigmask] = run_reliability(correctiontype,filename,data,ftbl);
        summary{curr_method,curr_GSR} = icc_summary{2,2}(end,end);
        summary_all{curr_method,curr_GSR} = icc_summary;

        save(['/home/kailong/Scheinost-Lab/icc_summary_' filename],'icc_summary','var','stats','sigmask')     
%         end
    end
end

save(['/home/kailong/cpm_reliability_summary_' filename],'summary','summary_all');
