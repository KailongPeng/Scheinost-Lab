function s = kailong_minus(a,b)
if contains(a,b)
    begin_index = strfind(a,b);
    end_index = begin_index + length(b) - 1;
    s = a;
    s(begin_index:end_index) = [];
end 
    