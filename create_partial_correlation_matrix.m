%CPM: connectome based predictive modeling
%NBS: network based statistics
clear all;close all; clc;
restoredefaultpath;
if isdir('/home/kailong/Scheinost-Lab')
    workdingdir = '/home/kailong/Desktop/';
    addpath('/home/kailong/Scheinost-Lab');
else
    workingdir = '/Users/pengkailong/Desktop/0 Yale/courses/rotation/Dustin Scheinost/';
    addpath('/Users/pengkailong/Desktop/0 Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab')
end

folder = [workdingdir 'results_matrix_268_110817'];
fileList = dir([folder '/*txt']);

%TRT005_2_TB_S006_bis_matrix_roimean.txt.txt
%TRT[subject]_[session]_[scanner]_S00[run]*
% t = extractfield(fileList,'name');
t = cell(1,size(fileList,1));
for curr_file = 1:size(fileList,1)
    t{curr_file} = fileList(curr_file).name;
end
fileList = [];
fileList = t;
t = [];
if ~isdir([folder '/partial_correlation'])
    mkdir([folder '/partial_correlation']);
end
parfor curr_file = 1:size(fileList,2)
    filename = [folder '/' fileList{curr_file}];
    savefile = [folder '/partial_correlation/' fileList{curr_file}(1:end-4)];
    if exist([savefile '.mat'])
%         movefile([savefile '.mat'],[newsavefile '.mat']);
        fprintf('file exist\n')
        continue;
    end
    curr_run_data = [];
    curr_run_data = readtable(filename);
    curr_run_data(:,end) = [];
    curr_run_data(:,1) = [];
    curr_run_data = table2array(curr_run_data);
    
    tic
    rho = partialcorr(curr_run_data);
    toc
    
    save(savefile,'rho');
    rho = [];
    filename = [];
    savefile = [];
end







