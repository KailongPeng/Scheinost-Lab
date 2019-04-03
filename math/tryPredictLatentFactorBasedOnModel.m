% X = zscore(X);
clear all;close all;
cd('/Users/pengkailong/Documents/MATLAB/Examples/stats/EstimateAndPlotFactorLoadingsExample')
m = 2;
load('X')

X = X(randperm(size(X,1)),:);
[Lambda,Psi,T,stats,F] = factoran(X,m,'scores','regression','rotate','none');%
X = X(1:100,:);
Fa = estimateFSBasedOnModel(X,Lambda,Psi,T,m);
Fa - F(100,:)
figure;plot(Fa(:,1),F(1:100,1),'.')
mdl = fitlm(Fa(:,1),F(1:100,1));
figure;boxplot(Fa)
figure;boxplot(F)
figure;boxplot(Fa - F(100,:))

function Fa = estimateFSBasedOnModel(X,Lambda,Psi,T,m,F)
% F can be empty
% example model:
% [Lambda,Psi,T,stats,F] = factoran(X,m,'scores','regression');

n = size(X,1);
stdev = std(X);
sqrtPsi = sqrt(Psi);
invsqrtPsi = diag(1 ./ sqrtPsi);
X0 = (X - repmat(mean(X),n,1)) ./ repmat(stdev,n,1);
Fa = [X0*invsqrtPsi zeros(n,m)] / [Lambda'*invsqrtPsi T'];
if exist('F')
    if ~isempty(F)
        if sum(sum(Fa - F)) ~= 0
            error('calculate F wrong!\n');
        end
    end
end
end
%
% L = Lambda;
% size(L)
% size(X)
% factor = inv(L'*L)*L'*(X'-mean(mean(X)));
% s = factor' - F ;
%
%
% F = F';
% u = factor(1,:);
% v = factor(2,:);
% CosTheta1 = dot(u,v)/(norm(u)*norm(v));
% u = F(1,:);
% v = F(2,:);
% CosTheta2 = dot(u,v)/(norm(u)*norm(v));
% u = factor(1,:);
% v = F(1,:);
% CosTheta3 = dot(u,v)/(norm(u)*norm(v));
%
% u = factor(2,:);
% v = F(2,:);
% CosTheta4 = dot(u,v)/(norm(u)*norm(v));
%


