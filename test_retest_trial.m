%{
%You should be able to find the originals here:
%/mnt/store1/mridata2/mri_group/smn33_data/test_retest/results/results_matrix_268_110817/

%Typically I would use "load_reliability_data" to load in the matrices
%and create factor tables from that, then give the data and factor 
%table to "run_reliability". Use "none" for the correction type.

clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab/'));
%pathname = '/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST2_LR/matrices/';
global kailongLog curr_log
kailongLog = [];
curr_log = 1;
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
GSR_list = {'GSR' 'NOGSR'};
% GSR_list = {'GSR'};

for curr_GSR = 1:size(GSR_list,2)
    GSR = GSR_list{curr_GSR};
    for curr_task = 1:size(task_list,2)
        clear task task_matrice
        task = task_list{curr_task};
        task_matrice = ['MxMxN_matrix_' task];
        savefolder = []; 
        savefolder = ['/home/kailong/Desktop/test_retest/' GSR '/'];
        if ~isdir(savefolder)
            mkdir(savefolder)
        end
        if exist([savefolder task_matrice '.mat'])
            fprintf('exist\n')
            continue
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
%         index = mat2cell(index);
%         index = index{~cellfun((@(x) isnan(x)),index)};
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
    end
end

for curr_log = 1:size(kailongLog,2)
    fprintf(num2str(kailongLog{curr_log}));
end
format long g
save(['/home/kailong/Scheinost-Lab/Log/' datestr(datetime) '_Log'],'kailongLog');
%}
%%
clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab/'));
% cd('/home/kailong/Scheinost-Lab/ICC_toolbox')
% show = 'show';
% save('/home/kailong/Scheinost-Lab/ICC_toolbox/show','show');
thisFolder = '/mnt/store1/mridata2/mri_group/smn33_data/test_retest/results/results_matrix_268_110817';
% thisFolder = '/home/kailong/Desktop/results_matrix_268_110817/';
cd(thisFolder);


prompt = 'What is the factors wanted? e.g. 1 2 3 ; 1 4 5 ; 4 5\n';
factors = input(char(prompt),'s');

% filename = 'test_retest_original_full_correlation';
% filename = 'test_retest_new_full_correlation';
% filename = 'test_retest_new_partial_correlation';
% filename = 'test_retest_new_full_correlation_from_dustin';
% filename = 'test_retest_new_full_correlation_no_atanh_from_dustin';
% filename = 'test_retest_new_full_correlation_no_atanh';
filename = 'test_retest_new_partial_correlation_no_atanh';
% partial = input('partial correlation? y or n\n','s');
% if strcmp(partial,'y')
if contains(filename,'partial')
    partial = 'partial_correlation';
else
    partial = '';
end

filename = [filename '_' strrep(factors,' ','_')];

if contains(filename,'original')
    thisPattern = '.*matrix_matrix\.txt$';
else
    thisPattern = '.*roimean\.txt';
end

if contains(filename,'dustin')
    thisFolder = '/home/kailong/Desktop/results_matrix_268_110817';
end


if ~exist(['/home/kailong/Desktop/' filename partial '.mat'])
    if strcmp(partial,'partial_correlation')
        [data,ftbl] = load_reliability_data_partial_correlation(thisFolder, thisPattern,filename);
    else
        [data,ftbl] = load_reliability_data(thisFolder, thisPattern,filename);
    end

    %ftbl: 1:subj(1-12) 2:scanner(1-2) 3:corrected session with same scanner 4:run(1-6) 5:session(1-4)
    prompt = 'What is the run number each time? ';
    % run_number = input(prompt);
    run_number = 6;
    ftbl(:,4) = repmat([1:run_number]',[size(ftbl,1)/run_number,1]);
%     save(['/home/kailong/Desktop/test_retest_trial_' partial],'data','ftbl');
    save(['/home/kailong/Desktop/' filename partial],'data','ftbl');
end
%%
% partial = input('partial correlation? y or n\n','s');
% if strcmp(partial,'y')
%     partial = 'partial_correlation';
% else
%     partial = '';
% end
% load(['/home/kailong/Desktop/test_retest_trial'],'data','ftbl');
% load(['/home/kailong/Desktop/test_retest_trial_' partial '.mat'],'data','ftbl');
load(['/home/kailong/Desktop/' filename partial '.mat'],'data','ftbl');

eval(['ftbl = ftbl(:,[' factors ']);']);
% ftbl(:,[4,5]) = [];
correctiontype = 'none';                                                                  
[icc_summary,var,stats,sigmask] = run_reliability(correctiontype,filename,data,ftbl);

save(['/home/kailong/Scheinost-Lab/icc_summary_' filename partial],'icc_summary','var','stats','sigmask')
%{
%check data format: find the runs that is not 6 min, e.g. shorter than 6
%min= 360 datapoints
matsize = cell2mat(cellfun((@(x) size(x)),data,'UniformOutput',0));
ID1 = find(matsize(1:2:end) ~= mode(matsize(1:2:end)));
ID2 = find(matsize(2:2:end) ~= mode(matsize(2:2:end)));
small_ID_to_delete = [ID1 ID2];

% for curr_ID2delete = size(small_ID_to_delete,2):-1:1
%     data{small_ID_to_delete(curr_ID2delete)}=[];
% end
% ftbl = ftbl(~cellfun('isempty',data),:);
% data = data(~cellfun('isempty',data));

% ftbl(find(ftbl(:,4)==7),1)
suj_to_delete = ftbl(small_ID_to_delete,1);
for ii = 1:size(suj_to_delete,1)
    ftbl(ftbl(:,1) == suj_to_delete(ii),:) = nan;
end
data = data(~isnan(ftbl(:,1)));
ftbl = ftbl(~isnan(ftbl(:,1)),:);
%}
