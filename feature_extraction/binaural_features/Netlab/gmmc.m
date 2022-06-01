function mixc = gmmc(nin, nClasses, ncentres, covarType)
% GMMC Creates GMM classification model

mixc.type = 'gmmc';
% Number of gaussians per class
mixc.nin  = nin;

% Loop over number of classes
for m = 1:nClasses
    mixc.model{m} = gmm(nin, ncentres, covarType);
end

% Initialize priors
mixc.priors = ones(1, nClasses)./nClasses;