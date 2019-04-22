%% preconfiguration

% address for LR matrix
path = '/data_dustin/math2/';

% all the dataset names
dataset = cell(4,1);
dataset{1} = 'Mult';
dataset{2} = 'Sub';
dataset{3} = 'Num';
dataset{4} = 'Rhyming';
num_task = 4;
num_run = 7;

for curr_ses = 1:2
    %% read data from sub_comm
    sub_list = 1 : 132;
    motion_para = zeros(length(sub_list), num_run);
    for i_sub = sub_list
        for j_task  = 1 : 3
            % run1
            file_name = sprintf([path, 'sub-%03d/ses-T1/func/sub-%03d_realign/REALIGN_sub-%03d_ses-T1_task-%s_run-01_bold_frametoframe.txt'], i_sub, i_sub, i_sub, dataset{j_task});
            if exist(file_name, 'file') == 2
                fileID = fopen(file_name,'r');
                a = fscanf(fileID,'%s %f %f');
                motion_para(i_sub, (j_task-1)*2+1) = a(end);
            end
            
            % run2
            file_name = sprintf([path, 'sub-%03d/ses-T1/func/sub-%03d_realign/REALIGN_sub-%03d_ses-T1_task-%s_run-02_bold_frametoframe.txt'], i_sub, i_sub, i_sub, dataset{j_task});
            if exist(file_name, 'file') == 2
                fileID = fopen(file_name,'r');
                a = fscanf(fileID,'%s %f %f');
                motion_para(i_sub, (j_task-1)*2+2) = a(end);
            end
        end
        file_name = sprintf([path, 'sub-%03d/ses-T1/func/sub-%03d_realign/REALIGN_sub-%03d_ses-T1_task-Rhyming_bold_frametoframe.txt'], i_sub, i_sub, i_sub);
        if exist(file_name, 'file') == 2
            fileID = fopen(file_name,'r');
            a = fscanf(fileID,'%s %f %f');
            motion_para(i_sub, num_run) = a(end);
        end
    end
    
    %% exclude data by threshold
    mean_thres = 0.1;
    max_thres = 0.15;
    bad_id = [];
    for i_sub = sub_list
        temp_motion = motion_para(i_sub, :);
        temp_motion(temp_motion==0) = [];
        if (mean(temp_motion)>mean_thres) || (max(temp_motion)>max_thres)
            bad_id = [bad_id, i_sub];
        end
    end
    motion_para_corrected = motion_para;
    motion_para_corrected(bad_id, :)=[];
    save(['/home/kailong/Scheinost-Lab/math/data/TwoSession/motion_ses' num2str(curr_ses) '.mat'], 'motion_para', 'bad_id')
end