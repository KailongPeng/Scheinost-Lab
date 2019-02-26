% partial = [0.259086020436656;
% 0.174635834133107;
% 0.315968073394861;
% 0.292052624670823];
% 
% full=[0.137944883590487;
% 0.162245364965930;
% 0.0899442164146991;
% 0.153784460595690;];

0.124500451175548
0.168676484237794
0.0802926182942722
0.137148888162656

clear sum full_no_max
load('/home/kailong/full_no_max_p3_.mat');
for i=1:4
    full_no_max{i} = sum{[7+i],7};
end
full_no_max = cell2mat(full_no_max);


clear sum full
load('/home/kailong/full_correlation_p3_.mat');
for i=1:4
    full{i} = sum{[7+i],7};
end
full = cell2mat(full);



clear sum partial
load('/home/kailong/partial_correlation_p3_.mat');
for i=1:4
    partial{i} = sum{[7+i],7};
end
partial = cell2mat(partial);

[h,p] = ttest(full,partial);
figure;
plot([1,1,1,1],full_no_max,'ro');
hold on;
plot([2,2,2,2],full,'bo');
plot([3,3,3,3],partial,'go');
xlim([0,4])

