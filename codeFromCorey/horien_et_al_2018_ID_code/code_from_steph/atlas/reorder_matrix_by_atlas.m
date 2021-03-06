function reordered_mat = reorder_matrix_by_atlas(orig_data,category)
% orig mat is nxixi; typically a cell mat of connectivity matrix
% e.g., reordered_mat = reorder_matrix_by_atlas(orig_mat,'subnetwork')
% Don't use this natively to summarize ROIs; domatrixsummary_XXX to organize

if iscell(orig_data)
         matdim=size(orig_data{1});

else matdim=size(orig_data);
end

if matdim(1) ~= matdim(2)
    ismatrix=0;
    warning('Input is vector, not matrix. If should be matrix, re-structure first.')
else
    ismatrix=1;
end

matdim=matdim(1);
load(sprintf('map%0.0f_%s.mat',matdim,category));
% if matdim==278
%      load(sprintf('map%0.0f.mat',matdim)); % need: map268.mat or map278.mat
%      map_old=map; clearvars map; map.oldroi=map_old(:,2);


% reorder
if ismatrix
    if iscell(orig_data)
        for n=1:length(orig_data)
            reordered_mat{n}=orig_data{n}(map.oldroi,map.oldroi);
        end
    else
         reordered_mat=orig_data(map.oldroi,map.oldroi);
    end
    
else
    if iscell(orig_data)
        for n=1:length(orig_data)
            reordered_mat{n}=orig_data{n}(map.oldroi);
        end
    else
         reordered_mat=orig_data(map.oldroi);
    end
end


end

