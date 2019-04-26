function t = kailong_cell2mat(M)
t = [];
for ii = 1:length(M)
    t(:,:,ii) = M{ii};
end