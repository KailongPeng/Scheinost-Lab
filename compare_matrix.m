%CPM: connectome based predictive modeling
%NBS: network based statistics
%TRT005_2_TB_S006_bis_matrix_roimean.txt
%TRT[subject]_[session]_[scanner]_S00[run]*

clear all;close all; clc;
restoredefaultpath;
if isdir('/home/kailong/Scheinost-Lab')
    workdingdir = '/home/kailong/Desktop/';
    addpath(genpath('/home/kailong/Scheinost-Lab'));
else
    workingdir = '/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/';
    addpath('/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab')
end
thisFolder = '/mnt/store1/mridata2/mri_group/smn33_data/test_retest/results/results_matrix_268_110817';
thisPattern = '.*matrix_matrix\.txt$';
theseFiles = regexpdir(thisFolder, thisPattern);
% theseFiles = kailong_extractfield(theseFiles,'')

folder = [workdingdir 'results_matrix_268_110817'];

fileList = dir([folder '/*txt']);

t = kailong_extractfield(fileList,'name');
fileList = [];
fileList = t;
t = [];
if ~isdir([folder '/correlation'])
    mkdir([folder '/correlation']);
end
sum = [];
for curr_file = 1:size(fileList,2)
    filename = [folder '/' fileList{curr_file}];
    savefile = [folder '/correlation/' fileList{curr_file}(1:end-4)];
    if exist([savefile '.mat'])
%         movefile([savefile '.mat'],[newsavefile '.mat']);
        fprintf('file exist\n')
        continue;
    end
    clear curr_run_data ;
    curr_run_data = readtable(filename);
    curr_run_data(:,end) = [];
    curr_run_data(:,1) = [];
    curr_run_data = table2array(curr_run_data);
    [cr,lgs] = xcorr(curr_run_data,10,'coeff');
    cr_max = max(abs(cr));
    [shift,corrected_cr_max] = findshift(cr,cr_max,lgs);
    cr_max = [];
    cr_max = corrected_cr_max;
    corrected_cr_max = [];
    cr_max = reshape(cr_max,[268,268]);  
    cr_no_max = xcorr(curr_run_data,0,'coeff');
    cr_no_max = reshape(cr_no_max,[268 268]);
    atanh_cr_no_max = atanh(cr_no_max-diag(diag(cr_no_max)));
    
    clear curr_run_data_original ;
    curr_run_data_original = readtable(theseFiles{curr_file});
    curr_run_data_original(:,end) = [];
%     curr_run_data_original(:,1) = [];
    curr_run_data_original = table2array(curr_run_data_original);
    curr_run_data_original = curr_run_data_original - diag(diag(curr_run_data_original));
    
    
    sum(curr_file,1) = curr_file    
    sum(curr_file,2) = kailong_isequal(cr_max,curr_run_data_original)
    sum(curr_file,3) = kailong_isequal(cr_no_max,curr_run_data_original)
    sum(curr_file,4) = kailong_isequal(atanh_cr_no_max,curr_run_data_original)

    %     save(savefile,'cr','cr_max','shift');
    cr = [];
    cr_max = [];
    shift = [];
    filename = [];
    savefile = [];
end
save('/home/kailong/Scheinost-Lab/compare_matrix_summary','sum')
