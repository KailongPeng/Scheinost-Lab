function [cr_max,shift] = create_max_correlation_matrix_function(filename)
curr_data = readtable(filename);
curr_data(:,end) = [];
curr_data(:,1) = [];
curr_data = table2array(curr_data);
[cr,lgs] = xcorr(curr_data,10,'coeff');
cr_max = max(abs(cr));
if ~isempty(find(isnan(cr_max)==1))
    error('nan')
end
[shift,corrected_cr_max] = findshift(cr,cr_max,lgs);
cr_max = [];
cr_max = corrected_cr_max;
corrected_cr_max = [];
cr_max = reshape(cr_max,[268,268]);
cr_max = cr_max - diag(diag(cr_max));
% cr_max = atanh(cr_max);
