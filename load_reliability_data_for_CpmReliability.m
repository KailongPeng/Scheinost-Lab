function [data,ftbl] = load_reliability_data_for_CpmReliability(thisFolder, thisPattern,filename,task_list,sub_num)
% built to parse metadata from traveling_subs and test_retest
% supports regex
% e.g., [data,ftbl] = load_trt_files('/Users/stephanie/Documents/data/mnt/gsr_unism/','.*roi_.\.nii\.gz');
% out:
%   ftbl(subj,scanner,day)
%   data = cell matrix; each cell contains an image vector corresponding with a file
% remember! for trav seeds, use /more_results/pcc/*crop_resampled.nii.gz

flag_studytype=0;

if ~isfolder(thisFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', thisFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% filePattern = fullfile(thisFolder, thisPattern);
% filePattern = strcat(thisFolder, thisPattern);
% theseFiles = dir(filePattern);
if contains(thisPattern,'cpm_reliability')
    for i = 1:length(task_list)
        theseFiles{i} = [thisFolder 'MxMxN_matrix_' task_list{i} '.mat'];
    end
else 
    theseFiles = regexpdir(thisFolder, thisPattern);
end

if isempty(theseFiles)
    errorMessage = sprintf('Error: Nothing found matching\n%s', thisPattern);
    uiwait(warndlg(errorMessage));
    return;
end

if contains(filename,'original')
%     thisPattern = '.*matrix_matrix\.txt$';
%     str = 'generate new correlation matrix? y or n\n';
    new_correlaltion = 'n';

else
    new_correlaltion = 'y';
end


for k = 1:length(theseFiles)
%     baseFileName = theseFiles(k).name;
%     fullFileName = fullfile(thisFolder, baseFileName);
      fullFileName = theseFiles{k};

          % report percent completion
        if ~exist('amt_fin_old')
            amt_fin_old=0;
        end
        if (mod(round(k*100/length(theseFiles)),5)==0)
            amt_fin=round(k*100/length(theseFiles));
            if(~(amt_fin==amt_fin_old))
                amt_fin_old=amt_fin;
                fprintf('%0.0f%% finished\n',amt_fin);
            end
        end

%          fprintf(1, 'Now reading %s\n', fullFileName);
      
    [pathstr,name,ext] = fileparts(fullFileName);

    baseFileName=[name,ext];

    % check correct filetype and assign study_type
    if ~flag_studytype
        if(~isempty(strfind(baseFileName,'000_')))
            study_type='traveling_subs';
        elseif (~isempty(strfind(baseFileName,'TRT')))
            study_type='test_retest';
        elseif contains(baseFileName,'MxMxN_matrix_')
            study_type = 'cpm_reliability';
        else
            study_type='undefined';
        end
        if ~(strfind(ext,'.nii')) & ~(strfind(name,'.nii'))
            errorMessage = sprintf('Error: The file is not .nii or .nii.gz:\n%s', fullFileName);
            uiwait(warndlg(errorMessage));
            return;
        end
        flag_studytype=1;
    end
  

    % if(strfind(ext,'.txt'))
    % don't use this. scripts built to use output of nii.
    % see final_project for info on squeezing out singletons
    %       a=readtable(fullFileName,'delimiter', 'tab');
    %       a(:,3) = [];
    %       output{k} = a; 
    %{
    tic
    clear temp;
    temp = readtable(fullFileName,'ReadVariableNames',false);
    temp(:,1) = [];
    temp(:,end) = [];
    temp(1,:) = [];
    temp = table2array(temp);
    temp = cellfun(@(x) str2double(x),temp);
    data{k} = temp;
    clear temp;
    toc
    %}
    % file prototype: MxMxN_matrix_REST2_LR.nii
    % ftbl: direction session
    identifier=strfind(baseFileName,'REST');
    %this_name=baseFileName((identifier+4):(identifier+15)); % returns "01_1_TA_S001"
    if contains(baseFileName,'LR')
        this_direction = '1';
    elseif contains(baseFileName,'RL')
        this_direction = '2';
    end
    if contains(baseFileName,'REST2')
        this_session = '2';
    elseif contains(baseFileName,'REST')
        this_session = '1';
    end

    % generate new correlation matrix
    if strcmp(new_correlaltion , 'y')
        clear cr_max;
        try
            [cr_max,~] = create_max_correlation_matrix_function(fullFileName);
            data{k} = cr_max;
        catch
            data{k} = nan;
        end
        clear cr_max;
    else
        %use generated matrix
        clear matrix_file
        if contains(filename,'cpm_reliability')
            eval(['clear MxMxN_matrix_' task_list{k} ';']);
            load(fullFileName,['MxMxN_matrix_' task_list{k}]);
            eval(['matrix_file = MxMxN_matrix_' task_list{k} ';']);
            for curr_sub = 1:sub_num
                data{(k-1)*sub_num+curr_sub} = matrix_file(:,:,curr_sub);
                ftbl((k-1)*sub_num+curr_sub,1) = curr_sub;
                ftbl((k-1)*sub_num+curr_sub,2)=str2double(this_direction);
                ftbl((k-1)*sub_num+curr_sub,3)=str2double(this_session); 
            end
            
        else
            matrix_file = readtable(fullFileName);
            matrix_file(:,end) = [];
            matrix_file = table2array(matrix_file);
            data{k} = matrix_file;
        end
    end
%     switch study_type
% %     if strcmp(study_type,'traveling_subs')
%         case'traveling_subs'
%             % file prototype: 01_V3000_2_bis_matrix.nii
%             % ftbl: subj scanner day (within a particular scanner)
%             identifier=strfind(baseFileName,'000_');
%             this_name=baseFileName((identifier-5):(identifier+4)); % returns "01_V1000_1"
%             this_subj=this_name(5);
%             this_scanner=this_name(2);
%             this_day=this_name(10);
% 
%             ftbl(k,1)=str2num(this_subj);
%             ftbl(k,2)=str2num(this_scanner);
%             ftbl(k,3)=str2num(this_day);
% 
% %     elseif strcmp(study_type,'test_retest')
%         case 'test_retest'
%             % file prototype: TRT001_1_TA_S001_bis_matrix.nii
%             % ftbl: subj scanner day (within a particular scanner) run session
%             identifier=strfind(baseFileName,'TRT');
%             this_name=baseFileName((identifier+4):(identifier+15)); % returns "01_1_TA_S001"
%             this_subj=this_name(1:2);
%             this_session=this_name(4);
%             this_run=this_name(12);
% 
%             % coding scanner 1=TA, 2=TB
%             if(strcmp(this_name(6:7),'TA')); this_scanner='1';
%             else this_scanner='2';
%             end
% 
%             this_day=this_name(10);
% 
%             ftbl(k,1)=str2num(this_subj);
%             ftbl(k,2)=str2num(this_scanner);
%             % day/occasion is added after reading all
%             if ~(strcmp(this_run,'_'))
%                 ftbl(k,4)=str2num(this_run);
%                 else ftbl(k,4)=1;
%             end
%             ftbl(k,5)=str2num(this_session);
%         case 'cpm_reliability'
            %{
            this_subj=this_name(1:2);
            this_session=this_name(4);
            this_run=this_name(12);

            % coding scanner 1=TA, 2=TB
            if(strcmp(this_name(6:7),'TA')); this_scanner='1';
            else this_scanner='2';
            end

            this_day=this_name(10);

            % day/occasion is added after reading all
            if ~(strcmp(this_run,'_'))
                ftbl(k,4)=str2num(this_run);
                else ftbl(k,4)=1;
            end
            ftbl(k,5)=str2num(this_session);

            %}
%         otherwise
%             ftbl=[];

    end
    
end

% assign day/occasion (corrected session with same scanner)
% if strcmp(study_type,'test_retest')
%     for subj=unique(ftbl(:,1))'
%         for scanner=unique(ftbl(:,2))'
%             
%             t=(find(ftbl(:,1)==subj & ftbl(:,2)==scanner));
%             t(:,2)=ftbl(t,5);
%             
%             l=unique(t(:,2));
%             l=sortrows(l);
%             for i=1:length(l)
%                 t(t(:,2)==l(i),2)=i;
%             end
%             ftbl(t(:,1),3)=t(:,2);
%             
%         end
%     end
% end
  

