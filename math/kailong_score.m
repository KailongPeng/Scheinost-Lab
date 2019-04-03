function score = kailong_score(data,coeff,mu)

% [coeff,score,latent,tsquared,explained,mu] = pca(data,'algorithm','als');
% t = score1*coeff' + repmat(mu1,13,1);


score = (data - repmat(mu,size(data,1),1))/coeff';



