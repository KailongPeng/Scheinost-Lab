function [y_predict,summary_feature]=cpm_test(x,mdl,pmask)
% Test a Connectome-based Predictive Model using previously trained model
% x            Predictor variable
% mdl          Coefficient fits for linear model relating summary features to y
% pmask        Mask for significant features
% y_predict    Predicted y values

% For each subject, create summary feature and use model to predict y
if strcmp(class(mdl),'double')
    for i=1:size(x,2)
        summary_feature(i)=nansum([nanmean(x(pmask>0,i)),-nanmean(x(pmask<0,i))]);
        y_predict(i)=mdl(2)*summary_feature(i) + mdl(1); 
    end
else
    if iscell(mdl)
        for i=1:size(x,2)
            summary_feature(i)=nansum([nanmean(x(pmask>0,i)),-nanmean(x(pmask<0,i))]);
            func = mdl{1};
            beta = mdl{2};
            y_predict(i)=func(beta,summary_feature(i)); 
        end
    end
end    