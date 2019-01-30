function [masked_data]=get_masked_data(data,mask,correctiontype)
% data is cell array 1xnsessions; mask is just one matrix

mask=logical(mask); % NEW

if strcmp(correctiontype,'none')
    data_new=[];
    for(i=1:length(data))
        data_new=data{i}(:);
        masked_data{i}=data_new;
    end
else
    data_new=[];
    for(i=1:length(data))
        data_new=data{i}(mask);
        masked_data{i}=data_new;
    end
end