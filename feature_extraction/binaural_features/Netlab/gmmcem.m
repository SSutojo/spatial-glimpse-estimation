function [mixc, options] = gmmcem(mixc, x, t, opt)
% GMMCEM  Train GMM classification model

% Transpose input 
x = transpose(x);

% Assume that data has 1-of-N encoding to set model priors
[ndata, nclasses] = size(t);
mixc.priors = sum(t)./nclasses;

for m = 1:nclasses
  % Extract data from mth class
  t_index = find(t(:, m) == 1);
  % Initialise mth GMM
  mixc.model{m} = gmminit(mixc.model{m}, x(t_index, :), opt{1});
  % Train mth GMM
  [mixc.model{m}, options] = gmmem(mixc.model{m}, x(t_index, :), opt{2});
end