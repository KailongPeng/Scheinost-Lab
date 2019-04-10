function EstimatedScore = estimateScoreBasedOnModel(y,coeff,mu,score1)
% score1 can be empty
% example model:
% [coeff1,score1,latent,tsquared,explained,mu] = pca(y,'algorithm','als');
EstimatedScore = (y - repmat(mu,size(y,1),1))/coeff';

if exist('score1','var')
    if ~isempty(score1)
%         figure;
%         plot(EstimatedScore(:,1),score1(:,1),'.');
%         [p,S,mu] = polyfit(EstimatedScore(:,1),score1(:,1),1);
%         y = polyval(p,EstimatedScore(:,1));
%         hold on; plot(EstimatedScore(:,1),y,'.r')
        if sum(sum(EstimatedScore - score1)) > 1e-6
            error('calculate Score wrong!\n');
        end
    end
end