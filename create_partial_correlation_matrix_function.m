function [rho,atanh_rho] = create_partial_correlation_matrix_function(filename)
curr_data = readtable(filename);
curr_data(:,end) = [];
curr_data(:,1) = [];
curr_data = table2array(curr_data);
%%
% use to detect missing nodes, maybe inefficient
[cr,~] = xcorr(curr_data,10,'coeff');
cr_max = max(abs(cr));
if ~isempty(find(isnan(cr_max)==1))
    error('nan')
end
%%
tic
rho = partialcorr(curr_data);
toc

rho = rho - diag(diag(rho));
atanh_rho = atanh(rho);
