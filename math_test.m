restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab'));
cd('/home/kailong/Scheinost-Lab')

%dustin data
%check for motion and the registration
cd('/data_dustin/math2/results/')
load('/data_dustin/math2/results/');

%siyuan data
cd('/mnt/store4/mri_group/siyuan_data/math');
load('/mnt/store4/mri_group/siyuan_data/math');

%have a look at the behavioral measures in /data_dustin/math2/phenotype/ 
%and see what looks the most interesting (the math scores may be a good 
%starting point). Once you have a plan, touch based with Siyuan/Javid for 
%the current mutli-dimensional CPM code and let’s see what we can predict 
%in this data set.
path = '/data_dustin/math2/phenotype/';
cd(path);
fnameList = dir([path '*.json']);
fnameList = kailong_extractfield(fnameList,'name');
clear task_description
for curr_task = 1:length(fnameList)
    fname = fnameList{curr_task};
    taskname = kailong_minus(fname,'.json');
    task_description{curr_task}.task = jsondecode(fileread(fname));%'adhd-rs.json';
    task_description{curr_task}.taskname = taskname;%'adhd-rs.json';
    fname = [];
    taskname = [];
end
filename = 'ses-T1/adhd-rs.tsv';
s = tdfread(filename);