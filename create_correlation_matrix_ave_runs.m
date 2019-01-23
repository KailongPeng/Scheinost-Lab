clear all;close all; clc;
restoredefaultpath;
addpath('/Users/pengkailong/Desktop/0 Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab');

folder = '/Users/pengkailong/Desktop/0 Yale/courses/rotation/Dustin Scheinost/results_matrix_268_110817';
kailongLog = [];
curr_log = 1;
% average across runs (one file for each session) 
if ~isdir([folder '/correlation_by_session'])
    mkdir([folder '/correlation_by_session']);
end
subjectList = [1:12];
sessionList = [1:4];
%check whether file number fit.
fileList = dir([folder '/*.txt']);
if size(fileList,1) ~= size(subjectList,2)*size(sessionList,2)*6
    error('file number not fit!')
end
for curr_subject = 1:size(subjectList,2)
    subject = subjectList(curr_subject);
    for curr_session = 1:size(sessionList,2)
        session = sessionList(curr_session);

        %TRT005_2_TB_S006_bis_matrix_roimean.txt
        subject_session = [];
        subject_session = dir([folder '/TRT' strrep(sprintf('%3d',subject),' ','0') '_' sprintf('%d',session)...
            '_T*_S*_bis_matrix_roimean.txt']);
            subject_session = extractfield(subject_session,'name');
            savefile = [folder '/correlation_by_session/' subject_session{1}(1:11)];
            if exist([savefile '.mat'])
                fprintf('file exist\n')
                continue;
            end
            curr_session_table = [];
            curr_session_mat = [];
            if size(subject_session,2)>6
                fprintf('more than 6 runs\n')
                kailongLog{curr_log} = ['more than 6 runs\n'];
                curr_log = curr_log + 1; 
            end
            for curr_run = 1:size(subject_session,2)
                filename = [];
                filename = [folder '/' subject_session{curr_run}];
                
                curr_session_table{curr_run} = readtable(filename);
                curr_session_table{curr_run}(:,end) = [];
                curr_session_table{curr_run}(:,1) = [];
                curr_session_table{curr_run} = table2array(curr_session_table{curr_run});
                if size(curr_session_table{curr_run},1) ~= 360
                    fprintf(['run duration is not 6 min! ' subject_session{curr_run} ' \n' ...
                        'this session contains ' num2str(size(subject_session,2)) ' runs\n'])
                    kailongLog{curr_log} = ['run duration is not 6 min! ' subject_session{curr_run} ' \n' ...
                        'this file contains ' num2str(size(subject_session,2)) ' runs\n'];
                    curr_log = curr_log + 1; 
                end
                try %if in some runs, time duration is not 6 min (360 rows), skip this run
                    curr_session_mat(:,:,curr_run) = curr_session_table{curr_run};
                end
            end
            if size(curr_session_table{1},1) ~= 360
                error('duration wrong');
            end

            curr_session_mean = nanmean(curr_session_mat,3);
            [cr,lgs] = xcorr(curr_session_mean,10,'coeff');

            cr_max = max(abs(cr));
            [shift,corrected_cr_max] = findshift(cr,cr_max,lgs);
            cr_max = [];
            cr_max = corrected_cr_max;
            corrected_cr_max = [];
            save(savefile,'cr','cr_max','shift');
            cr = [];
            cr_max = [];
            shift = [];
            filename = [];
            savefile = [];
    end
end
for curr_log = 1:size(kailongLog,2)
    fprintf(kailongLog{curr_log});
end
