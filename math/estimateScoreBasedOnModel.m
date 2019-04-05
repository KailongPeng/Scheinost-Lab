function score = estimateScoreBasedOnModel(y,coeff,mu,score1)
% score1 can be empty
% example model:
% [coeff1,score1,latent,tsquared,explained,mu] = pca(y,'algorithm','als');
EstimatedScore = (y - repmat(mu,size(y,1),1))/coeff';

if exist('score1','var')
    if ~isempty(score1)
        figure;
        plot(score(:,1),score1(:,1));
        
        if sum(sum(EstimatedScore - score1)) ~= 0
            error('calculate Score wrong!\n');
        end
    end
end