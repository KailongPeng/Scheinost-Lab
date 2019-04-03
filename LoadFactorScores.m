cd('/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab/R/math/math/output');
clear all;
DataName = 'noHighCor_norm_SelectedTest_';%noHighCor_norm_SelectedTest_1FactorScores
path = ['/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab/R/math/math/output/' ...
    DataName];
numOfFactorList = [1:7];

summary = [];
for curr_numOfFactor = 1:length(numOfFactorList)
    numOfFactor = numOfFactorList(curr_numOfFactor);
    text = [path num2str(numOfFactor) 'FactorScores.csv'];
    % [data, result] = readtext(text);
    % temp = char(data);
    data = [];
    data = readtable(text);
    data = data(:,2:end);
    data = table2array(data);
    save([path num2str(numOfFactor) 'LatentFactorEstimate'],'data')
end
