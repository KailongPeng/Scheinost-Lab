for ii=.001:.005:1
    temp = L1precisionBCD(cr,ii);
    temp = temp-diag(diag(temp));
    if ~isempty(find(temp~=0))
        pause
    end
end



for ii=0:.05:1
    tic    
    temp = temp-diag(diag(temp));
    temp = L1precisionBCD(cr,ii);
    toc
    if ~isempty(find(temp~=0))
        pause;
    end
end
old_temp = temp;


aaa=real(temp);