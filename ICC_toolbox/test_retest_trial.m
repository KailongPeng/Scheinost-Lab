%Typically I would use "load_reliability_data" to load in the matrices
%and create factor tables from that, then give the data and factor 
%table to "run_reliability". Use "none" for the correction type.

clear all;close all;clc;
restoredefaultpath;
addpath('/home/kailong/Scheinost-Lab/ICC_toolbox')

thisFolder = '/home/kailong/Desktop/results_matrix_268_110817/';
thisPattern = '.*roimean\.txt';
[data,ftbl] = load_reliability_data(thisFolder, thisPattern);
%ftbl: 1:subj(1-12) 2:scanner(1-2) 3:corrected session with same scanner 4:run(1-6) 5:session(1-4)
save(['/home/kailong/Desktop/test_retest_trial'],'data','ftbl');
correctiontype = 'none';
[icc_summary,var,stats,sigmask] = run_reliability(correctiontype,data,ftbl);

