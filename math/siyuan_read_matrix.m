% This code read all the data in the HCP dataset and store them in a n_node*n_node*n_sub*n_tasks cell
% Please configure the address before running the code.

%% preconfiguration
% clear all;
% address for LR matrix
path = '/data_dustin/math2/results/';


% all the dataset names
num_run = 7;

%% read data from sub_comm
sub_list = 1 : 132;
num_node = 268;
all_mats = zeros(num_node, num_node, length(sub_list), num_run);
incomp_id = [];
for i_sub = sub_list
    disp(i_sub)
    count = 0;
    for j_run  = 1 : num_run
        file_name = sprintf([path, 'sub-%03d/ses1/sub-%03d_bis_matrix_%d_matrix.txt'], i_sub, i_sub, j_run);
        if exist(file_name, 'file') == 2
            count = count + 1;
            all_mats(:, :, i_sub, j_run) = load(file_name);
        end
    end
    if count < 7
            incomp_id = [incomp_id, i_sub];
    end
end

save('/home/kailong/Scheinost-Lab/math/data/all_mats.mat', 'all_mats')
save('/home/kailong/Scheinost-Lab/math/data/incomp_id.mat', 'incomp_id')
