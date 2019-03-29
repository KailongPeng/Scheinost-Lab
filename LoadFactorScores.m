cd('/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab/R/math/math/output');
clear all;
%noHighCor_SelectedTest- 4 factor-efa-factor-loadings.csv
% text = '/Users/pengkailong/Downloads/Untitled.txt';
% path = ['/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab/R/math/math/output/' ...
%     'noHighCor_all_data- '];
% data = 'math';
DataName = 'noHighCor_SelectedTest_';%noHighCor_SelectedTest_1FactorScores
path = ['/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab/R/math/math/output/' ...
    DataName];
numOfFactorList = [1:7];

summary = [];
for curr_numOfFactor = 1:length(numOfFactorList)
    numOfFactor = numOfFactorList(curr_numOfFactor);
    text = [path num2str(numOfFactor) 'FactorScores.csv'];
    % [data, result] = readtext(text);
    % temp = char(data);
    
    data = readtable(text);
    data = data(:,2:end);
    data = table2array(data);
    save([path num2str(numOfFactor) 'LatentFactorEstimate'],mapID,data)
end
