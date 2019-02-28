function [r,p,pmask,mdl]=cpm_train(x,y,pthresh,corr_type,LinearFlag)
% Train a Connectome-based Predictive Model
% x            Predictor variable
% y            Outcome variable
% pthresh      p-value threshold for feature selection
% r            Correlations between all x and y
% p            p-value of correlations between x and y
% pmask        Mask for significant features
% mdl          Coefficient fits for linear model relating summary features to y

% Select significant features
[r,p]=corr(x',y,'type',corr_type);
pmask=(+(r>0))-(+(r<0)); 
pmask=pmask.*(+(p<pthresh));

% For each subject, summarize selected features 
for i=1:size(x,2)
    summary_feature(i)=nansum([nanmean(x(pmask>0,i)), - nanmean(x(pmask<0,i))]);
end

% Fit y to summary features
% mdl=robustfit(summary_feature,y','ols');  
if LinearFlag == 1
    if ~sum(~isnan(summary_feature))
        warning('pthresh too high');
        mdl=robustfit(summary_feature,y');
    else
        mdl=robustfit(summary_feature,y');
    end
%     figure;hold on;
%     plot(summary_feature,y','.')
%     x_temp = [min(summary_feature):(max(summary_feature)-min(summary_feature))/1000:max(summary_feature)];
%     plot(x_temp,mdl(2)*x_temp+mdl(1));
else
    [modelfun,beta] = kailong_nonlinear(summary_feature,y');
    mdl{1}=modelfun;
    mdl{2}=beta;
end

