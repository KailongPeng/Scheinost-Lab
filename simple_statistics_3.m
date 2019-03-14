addpath(genpath('/home/kailong/Scheinost-Lab'))
% path = '/home/kailong/Desktop/results_matrix_268_110817/cpm_nonlinear/GSR/full_correlation/';
path = '/home/kailong/Desktop/results_matrix_268_110817/cpm_nonlinear/NOGSR/partial_correlation/';
file_list = dir([ path '*.mat']);
file_list = kailong_extractfield(file_list,'name');
for curr_file = 1:length(file_list)
    file = file_list{curr_file};
    load([ path file])
    all_sum_r(:,:,curr_file) = sum_r;
end

[h,p] = ttest(all_sum_r(1,1,:),all_sum_r(1,2,:))