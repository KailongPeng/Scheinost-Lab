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
if exist('F','var')
    if ~isempty(F)
        % figure; plot(Fa(:,1),F(:,1),'.');
        if sum(sum(Fa - F)) ~= 0
            error('latent factor calculation wrong!\n');
        end
    end
end
end