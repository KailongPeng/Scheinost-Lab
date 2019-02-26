function [atanh_cr,cr] = create_no_max_correlation_matrix(filename)
curr_data = readtable(filename);
curr_data(:,end) = [];
curr_data(:,1) = [];
curr_data = table2array(curr_data);
[cr,~] = xcorr(curr_data,10,'coeff');
cr_max = max(abs(cr));
if ~isempty(find(isnan(cr_max)==1))
    error('nan')
end
cr = [];
[cr,~] = xcorr(curr_data,0,'coeff');

cr = reshape(cr,[268,268]);
atanh_cr = cr - diag(diag(cr));
cr = atanh(atanh_cr);
