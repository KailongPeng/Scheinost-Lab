function [ FList ] = deleteAllLockFiles(DataFolder) 


% this function find the files whose name contain extList and delete them.


findtext = '_lock.mat';
if nargin < 1 
    if isdir('/home/kailong/Desktop') == 1 
        DataFolder = '/home/kailong/Desktop/results_matrix_268_110817/partial_correlation';
%     else
%         DataFolder = '/gpfs/ysm/home/kp578/git/ML_spontaneous_activity/'; 
    end
end

DirContents=dir(DataFolder); 
FList=[];

if ~isunix 
NameSeperator='\'; 
else isunix 
NameSeperator='/'; 
end

 
% Here 'peg' is written for .jpeg and 'iff' is written for .tiff 
for i=1:numel(DirContents) 
    if(~(strcmpi(DirContents(i).name,'.') || strcmpi(DirContents(i).name,'..'))) 
        if(~DirContents(i).isdir) 
            extension=DirContents(i).name; 
            if(numel(find(strfind(extension,findtext)))~=0) 
                FList={[DataFolder,NameSeperator,DirContents(i).name]}; 
                delete(FList{1}); 
            end 
        else 
            getlist=deleteAllLockFiles([DataFolder,NameSeperator,DirContents(i).name]); 
            FList=cat(1,FList,getlist); 
        end 
    end
end 