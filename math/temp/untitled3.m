[Lambda,Psi,T,stats,F] = factoran(X,2,'scores','regression');
inv(T'*T);   % Estimated correlation matrix of F, == eye(2)
Lambda*Lambda' + diag(Psi); % Estimated correlation matrix
Lambda*inv(T);              % Unrotate the loadings
F*T';                       % Unrotate the factor scores


size(X) %392*5 X is an n-by-d matrix where each row is an observation of d variables. 
size(Lambda)  %5*2 maximum likelihood estimate, lambda, of the factor loadings matrix.
size(Psi) %5*1  maximum likelihood estimates of the specific variances as a column vector psi of length d. 
size(T) %2*2  the m-by-m factor loadings rotation matrix T.
size(F) % 392*2 predictions of the common factors, known as factor scores. F is an n-by-m matrix where each row is a prediction of m common factors.


loading = Lambda*inv(T);
factorScore = F*T';

close all;
figure;boxplot((F))

figure;boxplot((X*Lambda - F))
figure;boxplot((Lambda'*X')' - F)

figure;boxplot(((X + repmat(Psi',size(X,1),1))/Lambda'- F))


%%
close all;
figure;imagesc(loading)
figure;boxplot((normalize(X)))
figure;boxplot((F))
figure;boxplot((factorScore))

figure;boxplot((normalize(X)*loading - factorScore))
figure;boxplot((loading'*normalize(X)')' - factorScore)

figure;boxplot(((normalize(X) + repmat(Psi',size(X,1),1))/loading'- factorScore))
figure;boxplot(((normalize(X) - repmat(Psi',size(X,1),1))/loading')-normalize(F))



L = loading;
Y = X;
miu = Psi;
figure;boxplot((inv(L'*L)*L'*(Y-repmat(miu',size(X,1),1))')')
figure;boxplot((inv(L'*L)*L'*(Y-repmat(mean(Y),size(X,1),1))')')

fhat = (inv(L'*L)*L'*(Y-repmat(mean(Y),size(X,1),1))')';

fhat = (inv(L'*L)*L'*(Y-repmat(miu',size(X,1),1))')';


L = Lambda;

fhat = (inv(L'*L)*L'*(Y-repmat(mean(Y),size(X,1),1))')';

fhat = (inv(L'*L)*L'*(Y-repmat(miu',size(X,1),1))')';



figure;boxplot(fhat)




