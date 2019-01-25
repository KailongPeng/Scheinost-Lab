function a = kailong_extractfield(S,name)

for curr_struct = 1:size(S,1)
    eval(['a{curr_struct} = S(curr_struct).' name ';']);
end
