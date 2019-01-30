%You should be able to find the originals here:
%/mnt/store1/mridata2/mri_group/smn33_data/test_retest/results/results_matrix_268_110817/

%Typically I would use "load_reliability_data" to load in the matrices
%and create factor tables from that, then give the data and factor 
%table to "run_reliability". Use "none" for the correction type.

clear all;close all;clc;
restoredefaultpath;
addpath(genpath('/home/kailong/Scheinost-Lab/'));
% cd('/home/kailong/Scheinost-Lab/ICC_toolbox')
% show = 'show';
% save('/home/kailong/Scheinost-Lab/ICC_toolbox/show','show');
thisFolder = '/mnt/store1/mridata2/mri_group/smn33_data/test_retest/results/results_matrix_268_110817';
% thisFolder = '/home/kailong/Desktop/results_matrix_268_110817/';
cd(thisFolder);
thisPattern = '.*roimean\.txt';
[data,ftbl] = load_reliability_data(thisFolder, thisPattern);
%ftbl: 1:subj(1-12) 2:scanner(1-2) 3:corrected session with same scanner 4:run(1-6) 5:session(1-4)
prompt = 'What is the run number each time? ';
run_number = input(prompt);
ftbl(:,4) = repmat([1:run_number]',[size(ftbl,1)/run_number,1]);
save(['/home/kailong/Desktop/test_retest_trial_2'],'data','ftbl');

clear all;close all;clc;
% load(['/home/kailong/Desktop/test_retest_trial'],'data','ftbl');
load(['/home/kailong/Desktop/test_retest_trial_2.mat'],'data','ftbl');
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

ftbl(:,3) = [];

correctiontype = 'none';
[icc_summary,var,stats,sigmask] = run_reliability(correctiontype,data,ftbl);

