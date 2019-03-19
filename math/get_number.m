function num_id = get_number(str)
%     str
    for i = 1:size(str,1)
        if isempty(regexp(str(i,:),'-?+\d+(\.\d+)?','Match'))
            num_id(i) = nan;
        else
            try
                num_id(i) = str2double(regexp(str(i,:),'-?+\d+(\.\d+)?','Match'));
            catch
                if length(str2double(regexp(str(i,:),'-?+\d+(\.\d+)?','Match'))) > 1
%                     temp = [];
%                     temp = str2double(regexp(str(i,:),'-?+\d+(\.\d+)?','Match'));
                    num_id(i) = nan;
                end
            end
        end
    end
    if sum(isnan(num_id)) > size(str,1)/2
        num_id = nan([size(str,1),1]);
    end
end