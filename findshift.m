function [shift,corrected_cr_max] = findshift(cr,cr_max,lgs)

shift = nan([1,size(cr,2)]);
for curr_roi = 1:size(cr,2)
    row = [];
    [row,~]= find(cr(:,curr_roi) == cr_max(1,curr_roi));
    if isempty(row)
        cr_max(1,curr_roi) = - cr_max(1,curr_roi);
        [row,~]= find(cr(:,curr_roi) == cr_max(1,curr_roi));
    end
    shift(1,curr_roi) = lgs(row);
end
corrected_cr_max = [];
corrected_cr_max = cr_max;