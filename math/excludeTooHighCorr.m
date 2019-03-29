function [noHighCor_table,noHighCor_table_VariableNames,noHighCor_norm_table] = excludeTooHighCorr(cr,Variable_class,char_VariableNames,norm_no_nan_sum_all_test,no_norm_no_nan_sum_all_test)
% tooHighThreshList = [0.4 0.60 0.8 0.9];
tooHighThreshList = [0.8];
for curr_tooHighThresh = 1:length(tooHighThreshList)
    tooHighThresh = tooHighThreshList(curr_tooHighThresh);

    [row,col] = find(cr<tooHighThresh);
    cr_high = cr;
    for curr_point = 1:length(row)
        cr_high(row(curr_point),col(curr_point)) = nan;
    end
    cr_high = tril(cr_high);
    figure
    imagesc(cr_high)
    hold on
    diff_ID = diff(Variable_class);diff_ID = [diff_ID 1];
    t = find(diff_ID ~= 0);t = [0 t];
    l_list = diff(t);
    
%     figure
%     imagesc(cr)
%     hold on
% draw square
    for curr_squre = 1:length(l_list)
    rectangle('Position',[sum(l_list(1:curr_squre-1))+0.5,sum(l_list(1:curr_squre-1))+0.5,l_list(curr_squre),l_list(curr_squre)],...
              'Curvature',0,...
             'LineWidth',2,'LineStyle','-','EdgeColor','r')
    end
%     write test name
    for curr_variable = 1:size(char_VariableNames,1)
        if and(curr_variable < 30,size(char_VariableNames,1)>50)
            xt = 62;
        else
            xt = 1;
        end
        yt = curr_variable;
        str = strrep(char_VariableNames{curr_variable},'_',' ');
        text(xt,yt,str,'Color','red','FontSize',9)
    end  
    % delete the one of the pairs whose correlation is bigger than tooHighThresh to try
    % factor analysis again.
    [new_row,new_col] = find((and(~isnan(cr_high),cr_high~=0))==1);
    remain_ID = [1:size(norm_no_nan_sum_all_test,2)];
    remain_ID(unique(new_row))=[];
    noHighCor_norm_no_nan_sum_all_test = norm_no_nan_sum_all_test(:,remain_ID);
    noHighCor_no_norm_no_nan_sum_all_test = no_norm_no_nan_sum_all_test(:,remain_ID);
    char_VariableNames = char(char_VariableNames);
    char_VariableNames = char_VariableNames(remain_ID,:);
    char_VariableNames = cellstr(char_VariableNames);
    
%     FA_for_all_data(noHighCor_norm_no_nan_sum_all_test);
    noHighCor_norm_table = [];
    noHighCor_norm_table = array2table(noHighCor_norm_no_nan_sum_all_test,'VariableNames',char_VariableNames);
    
    noHighCor_table = [];
    noHighCor_table = array2table(noHighCor_no_norm_no_nan_sum_all_test,'VariableNames',char_VariableNames);
    

    
%     writetable(noHighCor_table,'/home/kailong/Scheinost-Lab/math/data/noHighCor_all_test_no_norm_no_nan_header','Delimiter',',')
    noHighCor_table_VariableNames = cell2table(char_VariableNames');
%     writetable(noHighCor_table_VariableNames,'/home/kailong/Scheinost-Lab/math/data/VariableNames','Delimiter',',')

end