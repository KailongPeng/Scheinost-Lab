function FA_for_all_data(noHighCor_norm_no_nan_sum_all_test)
% load('/home/kailong/Scheinost-Lab/math/data/all_test_no_norm_no_nan_header')
% norm_no_nan_sum_all_test = normalize(table2array(table),1);
% no_norm_no_nan_sum_all_test = table2array(table);


% %check the correlation between different tests
% [cr,~] = xcorr(norm_no_nan_sum_all_test,0,'coeff');
% cr = reshape(cr,[sqrt(length(cr)),sqrt(length(cr))]);
% cr = cr - diag(diag(cr));
% figure
% imagesc(cr)
% temp = tril(cr);
% temp(temp==0) = nan;
% figure;hist(temp(:),100);
% figure;boxplot(temp(:))
%%
%factor analysis
a = [];
for num_common_factor = 1:10 %when num_common_factor = 7, the warning Some unique variances are zero: cannot compute significance appears
    num_common_factor
    [a{num_common_factor}.lambda,a{num_common_factor}.psi,a{num_common_factor}.T,...
        a{num_common_factor}.stats,a{num_common_factor}.F] = ...
        factoran(noHighCor_norm_no_nan_sum_all_test,num_common_factor);
%     a{num_common_factor}.stats.p 
end
close all;
for num_common_factor = 1:10
    figure;
    imagesc(a{1, num_common_factor}.lambda>0.4)
    pause;
end
