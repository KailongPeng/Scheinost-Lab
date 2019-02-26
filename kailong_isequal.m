function equal = kailong_isequal(a,b)
c = a - b;
c(c<1e-05) = 0;
equal = isempty(find(c(:) ~= 0));