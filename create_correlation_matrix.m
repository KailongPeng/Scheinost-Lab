%CPM: connectome based predictive modeling
%NBS: network based statistics
%TRT005_2_TB_S006_bis_matrix_roimean.txt
%TRT[subject]_[session]_[scanner]_S00[run]*

clear all;close all; clc;
restoredefaultpath;
if isdir('/home/kailong/Scheinost-Lab')
    workdingdir = '/home/kailong/Desktop/';
    addpath('/home/kailong/Scheinost-Lab');
else
    workingdir = '/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/';
    addpath('/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab')
end

folder = [workdingdir 'results_matrix_268_110817'];

fileList = dir([folder '/*txt']);

t = extractfield(fileList,'name');
fileList = [];
fileList = t;
t = [];
if ~isdir([folder '/correlation'])
    mkdir([folder '/correlation']);
end
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
%     cr = reshape(cr,[268,268]);
    save(savefile,'cr','cr_max','shift');
    cr = [];
    cr_max = [];
    shift = [];
    filename = [];
    savefile = [];
end







