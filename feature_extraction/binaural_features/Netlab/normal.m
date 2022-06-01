function y = normal(x)
% NORMAL  Calculate normalisation to unit variance zero mean
%
% Linearly rescale data to have zero mean and unit variance
% in each column (apart from those with zero variance).
% 
% Copyright (c) Ian T Nabney (2000)

n = size(x, 1);
% Compute mean and variance
mu = mean(x, 1);
sigma = std(x, 1);

% Rescale data, taking care over columns with zero variance to avoid
% division by zero
e = ones(n, 1);
y = (x - e*mu);  % Make y have zero mean
y = y./(e*(sigma+(sigma==0)));


