% function compare_CPM_r_values

clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab/'));

% global kailongLog curr_log
% kailongLog = [];
% curr_log = 1;
task_list = {'REST_LR' 'REST_RL' 'REST2_LR' 'REST2_RL'};

GSR_list = {'GSR' 'NOGSR'};
% method_list = {'partial_correlation' 'full_correlation'};
method_list = {'partial_correlation' 'full_correlation' 'full_no_max'};
% method_list = {'full_correlation'};

for curr_method = 1:length(method_list)
    method = method_list{curr_method};
    for curr_GSR = 1:size(GSR_list,2)
        GSR = GSR_list{curr_GSR};
        for curr_task = 1:size(task_list,2)
            clear task task_matrice
            task = task_list{curr_task};
            task_matrice = ['MxMxN_matrix_' task];
            savefolder = []; 
    %         savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/partial_correlation/'];
    %         target_file = [savefolder task_matrice '_atanh.mat'];
            switch method
                case 'partial_correlation'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/partial_correlation/'];
                case 'full_correlation'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/'];
                case 'full_no_max'
                    savefolder = ['/home/kailong/Desktop/results_matrix_268_110817/' GSR '/full_no_max/'];
            end
            target_file = [savefolder task_matrice '_atanh_p3_.mat'];

            clear summary_r
            load([target_file],'summary_r','summary_p');
            
            sum{1,1} = 'task';
            sum{1,2} = 'max';
            sum{1,3} = 'varience';
            sum{1,4} = 'average';
            sum{1,5} = 'summary_r';
            sum{1,6} = 'summary_p';
            sum{1,7} = 'p0.05_k_max';
%             sum{1,7} = 'p_k';
            
            %task
            sum{curr_task+1+6*(curr_GSR-1),1} = task;
            %max
            [~,index] = nanmax(abs(summary_r(:)));
            sum{curr_task+1+6*(curr_GSR-1),2} = summary_r(index);
            %varience
            sum{curr_task+1+6*(curr_GSR-1),3} = nanvar(summary_r(:));
            %average
            sum{curr_task+1+6*(curr_GSR-1),4} = nanmean(summary_r(:));
            %original maxtrix
            sum{curr_task+1+6*(curr_GSR-1),5} = summary_r;
            %summary_p matrix
            sum{curr_task+1+6*(curr_GSR-1),6} = summary_p;
            
            sum{curr_task+1+6*(curr_GSR-1),7} = summary_r(2,3);
            
        end
    end
    save(['/home/kailong/' method '_p3_'],'sum')
    clear sum
end