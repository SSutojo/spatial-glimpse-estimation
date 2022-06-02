function [y, p] = gmmcfwd(mixc, x)
% GMMCFWD  Run GMM classifier
x = transpose(x);

ndata = size(x, 1);
nclasses = size(mixc.priors, 2);
y = zeros(ndata, nclasses);
if nargout == 2
  p = zeros(size(y));
end

% Run forward each model in turn
for m = 1:nclasses
    y(:, m) = gmmprob(mixc.model{m}, x)*mixc.priors(m);
end

if nargout == 2
  p = y./(sum(y, 2)*ones(1, nclasses));
end
