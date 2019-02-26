clear all;close all;clc;

load('/home/kailong/cpm_reliability_summary.mat')

for curr_method = 1:3
    for curr_GSR = 1:2
        temp{curr_method,curr_GSR} = summary_all{curr_method,curr_GSR}{2,2};
    end
end

% temp_all = temp{4};

for ii = 1:6
    for jj = 1:6
        [t_result_h(ii,jj),t_result_p(ii,jj)] = ttest(temp{ii}(:),temp{jj}(:));
    end
end

for ii = 1:6
    for jj = 1:6
%         [t_result_h(ii,jj),t_result_p(ii,jj)] = ttest([temp{ii}(:) - temp{jj}(:)],0,'Alpha',0.01);
        [t_result_h(ii,jj),t_result_p(ii,jj)] = ttest([temp{ii}(:) - temp{jj}(:)],0);
    end
end
for ii = 1:6
    for jj = 1:6
        minus_result(ii,jj) = mean(temp{ii}(:)) - mean(temp{jj}(:));
    end
end

minus_result(minus_result>0) = 1;
minus_result(minus_result<0) = 0;

