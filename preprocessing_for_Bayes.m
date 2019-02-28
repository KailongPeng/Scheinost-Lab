fileList = dir('/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/*.txt');
fileList = kailong_extractfield(fileList,'name');
filepath = '/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/';
for curr_file = 1:length(fileList)
    filename = [filepath fileList{curr_file}];
    writingname = [kailong_minus(filename,'.txt') '_1.txt'];
    curr_data = readtable(filename);
    curr_data(:,end) = [];
    curr_data(:,1) = [];
% curr_data = table2array(curr_data);
    writetable(curr_data,writingname,'Delimiter',' ')  
% xlswrite(writingname,curr_data)
% csvwrite(writingname,curr_data)
end
