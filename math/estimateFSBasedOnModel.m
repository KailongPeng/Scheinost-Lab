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