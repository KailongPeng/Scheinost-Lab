function [y_predict_reshape,summary_feature_reshape]=cpm_cv(x,y,pthresh,kfolds,corr_type,LinearFlag)
% Runs cross validation for CPM
% x            Predictor variable
% y            Outcome variable
% pthresh      p-value threshold for feature selection
% kfolds       Number of partitions for dividing the sample
% y_test       y data used for testing
% y_predict    Predictions of y data used for testing

% Split data
nsubs=size(x,2);
nfeats=size(x,1);
randinds=randperm(nsubs);
ksample=floor(nsubs/kfolds);

% Run CPM over all folds
fprintf('\n# Running over %1.0f Folds.\nPerforming fold no. ',kfolds);
y_test = [];
y_predict = [];
summary_feature = [];
for leftout = 1:kfolds
    fprintf('%1.0f ',leftout);
    
    if kfolds == nsubs % doing leave-one-out
        testinds=randinds(leftout);
        traininds=setdiff(randinds,testinds);
    else
        si=1+((leftout-1)*ksample);
        if leftout == kfolds
            fi=nsubs;
        else
            fi=si+ksample-1;
        end
        
        testinds=randinds(si:fi);
        traininds=setdiff(randinds,testinds);
    end
    nsubs_in_fold=length(testinds);
    
    % Assign x and y data to train and test groups 
    x_train = x(:,traininds);
    y_train = y(traininds);
    x_test = x(:,testinds);
%     y_test(leftout,1:nsubs_in_fold) = y(testinds);
    y_test = [y_test ; y(testinds)];
    
    % Train Connectome-based Predictive Model
    [r,p,pmask,mdl] = cpm_train(x_train, y_train,pthresh,corr_type,LinearFlag);
    
    % Test Connectome-based Predictive Model
%     [y_predict(leftout,1:nsubs_in_fold)]=cpm_test(x_test,mdl,pmask);
    temp = [];
    temp_summary_feature = [];
    [temp,temp_summary_feature]=cpm_test(x_test,mdl,pmask);
    y_predict = [y_predict temp];
    summary_feature = [summary_feature temp_summary_feature];
%     hist(temp_summary_feature);
end

% y_test_reshape(randinds)=reshape(y_test',[],1);
% y_predict_reshape(randinds)=reshape(y_predict',[],1);
y_test_reshape(randinds) = y_test;
y_predict_reshape(randinds) = y_predict';
summary_feature_reshape(randinds) = summary_feature';

