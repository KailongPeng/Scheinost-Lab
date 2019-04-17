function a = kailong_extractfield(S,name)
SClass = class(S);
switch SClass
    case 'cell'
        for curr_struct = 1:size(S,2)
            eval(['a{curr_struct} = S{curr_struct}.' name ';']);
        end
        
    case 'struct'
        for curr_struct = 1:size(S,1)
            eval(['a{curr_struct} = S(curr_struct).' name ';']);
        end
        
end