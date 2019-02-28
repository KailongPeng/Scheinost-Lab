corr_type_List = {'Pearson' 'spearman'};
LinearFlagList = [1,0];
curr_loop = 1;
p = 0.05;
% k_fold = size(MxMxN_matrix_temp,3);
kfold = 10;
for curr_corr_type = 1:length(corr_type_List)
    corr_type = corr_type_List{curr_corr_type};
    for curr_LinearFlag = 1:length(LinearFlagList)
        LinearFlag = LinearFlagList(curr_LinearFlag);
        [y{curr_loop}.y_predict, y{curr_loop}.performance] = cpm_main(MxMxN_matrix_temp,behavior(map_ID),corr_type,LinearFlag,'pthresh',p,'kfolds',kfold);
        sum_r(curr_corr_type,curr_LinearFlag) = y{curr_loop}.performance(1);
        sum_p(curr_corr_type,curr_LinearFlag) = y{curr_loop}.performance(2);
        curr_loop = curr_loop + 1;
    end
end
sum_r