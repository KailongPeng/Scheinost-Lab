sum_all_test = [];
clear VariableNames;
curr_VariableNames = 1;
for curr_task = 1:length(task_data)
    curr_task_fields_list = fields(task_data(curr_task).sesT1);
    for curr_task_field = 1:length(curr_task_fields_list)
        task_field = curr_task_fields_list{curr_task_field};
        if ~contains(task_field,'participant_id')
            if strcmp(class(eval(['task_data(curr_task).sesT1.' task_field ''])),'char')
                temp = [];
                temp = eval(['get_number(task_data(curr_task).sesT1.' task_field ')']);
                sum_all_test = [sum_all_test temp(:)];
            else
                temp = [];
                temp = eval(['task_data(curr_task).sesT1.' task_field ';']);
                sum_all_test = [sum_all_test temp(:)];
            end
            VariableNames{curr_VariableNames} = task_field;
            curr_VariableNames = curr_VariableNames + 1;
        end
    end
end
char_VariableNames = char(VariableNames);
char_VariableNames = char_VariableNames(logical(nansum(sum_all_test)),:);
char_VariableNames = cellstr(char_VariableNames);
no_nan_sum_all_test = sum_all_test(:,logical(nansum(sum_all_test)));
no_nan_sum_all_test = no_nan_sum_all_test(all(~isnan(no_nan_sum_all_test),2),:);
% sum(sum(isnan(no_nan_sum_all_test)))
norm_no_nan_sum_all_test = normalize(no_nan_sum_all_test,1);
table = array2table(norm_no_nan_sum_all_test,'VariableNames',char_VariableNames);
save('/home/kailong/Scheinost-Lab/math/data/all_test_norm_no_nan','table')
writetable(table,'/home/kailong/Scheinost-Lab/math/data/all_test_norm_no_nan','Delimiter',',')
% csvwrite('/home/kailong/Scheinost-Lab/math/data/all_test_norm_no_nan',table)

