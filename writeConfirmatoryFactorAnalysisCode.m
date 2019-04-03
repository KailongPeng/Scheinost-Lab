cd('/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab/R/math/math/output');
clear all;close all;

% DataName = 'noHighCor_SelectedTest';
DataName = 'noHighCor_norm_SelectedTest';
path = ['/Users/pengkailong/Desktop/Yale/courses/rotation/Dustin Scheinost/Scheinost-Lab/R/math/math/output/' ...
    DataName '- '];
numOfFactorList = [1:7];
% LoadingThreshold = 0.3;
LoadingThreshold = 0.5;
% LoadingThreshold = -10^9;
summary = [];
for curr_numOfFactor = 1:length(numOfFactorList)
    numOfFactor = numOfFactorList(curr_numOfFactor);
    text = [path num2str(numOfFactor) ' factor-efa-factor-loadings.csv'];
    data = readtable(text);
    VariableNames = data(:,1);
    VariableNames = table2cell(VariableNames);
    for curr_row = 1:size(VariableNames,1)
        temp{curr_row} = strrep(VariableNames{curr_row},' ','');
    end
    
    VariableNames = char(temp);
    VariableNames
    matrix = table2array(data(:,2:end));
    factor_component = matrix>LoadingThreshold;
    figure;imagesc(factor_component);
    
    all_str = [];
    for curr_factor = 1:size(factor_component,2)
        str = ['factor' num2str(curr_factor) ' =~ '];
        for ii = 1:sum(factor_component(:,curr_factor))
            list = find(factor_component(:,curr_factor)==1);
            str = [str temp{list(ii)} ' + '];
        end
        str(end-2:end) = [];
        all_str{curr_factor} = str;
        str = [];
    end
    summary{curr_numOfFactor} = all_str;
end
fileID = fopen('matlabWrittenRCode.txt','w');
for curr_numOfFactor = 1:length(numOfFactorList)
    numOfFactor = numOfFactorList(curr_numOfFactor);
    fprintf(fileID,['# ' num2str(numOfFactor) ' factor model \nmodels$m' num2str(numOfFactor) ' <- \n    '''])
    if numOfFactor ~= 1
        fprintf(fileID,['' summary{numOfFactor}{1} '\n\n'])
    end
    for curr_factor = 2:(numOfFactor-1)
        fprintf(fileID,['    ' summary{numOfFactor}{curr_factor} '\n\n'])
    end
    fprintf(fileID,['    ' summary{numOfFactor}{numOfFactor} '''\n\n\n\n'])
end
for curr_numOfFactor = 1:length(numOfFactorList)
    numOfFactor = numOfFactorList(curr_numOfFactor);
    fprintf(fileID,['fits$m' num2str(numOfFactor) ' <- lavaan::cfa(models$m' num2str(numOfFactor) ', data = ' DataName ')\n'])
end
% for curr_numOfFactor = 1:length(numOfFactorList)
%     numOfFactor = numOfFactorList(curr_numOfFactor);
%     fprintf(fileID,['fits$m' num2str(numOfFactor) ' <- lavaan::cfa(models$m' num2str(numOfFactor) ', data = ' DataName ')\n'])
% end
fclose(fileID);
% table2array(data)
% temp = table2cell(data);

% a = cellfun(( @(x) isempty(x)),data);
% index = find(a==1);
%
%
% hold on;
% % Then, from the help:
% figure
% rectangle('Position',[0.59,0.35,3.75,1.37],...
%           'Curvature',0.2,...
%          'LineWidth',2,'LineStyle','-','EdgeColor','r')
%
%      index = [1,6;7,14;15,23;24,33;34,69;70,75;76,77;78,82;83,96;97,107];
%

%                 Factor1 Factor2 Factor3 Factor4 Factor5 Factor6 Factor7 Factor8 Factor9
% SS loadings      6.828   3.719   3.663   3.607   3.544   3.078   3.003   2.596   2.184
% Proportion Var   0.100   0.055   0.054   0.053   0.052   0.045   0.044   0.038   0.032
% Cumulative Var   0.100   0.155   0.209   0.262   0.314   0.359   0.404   0.442   0.474