function a = loggmmactiv(mix, x)
% function a = loggmmactiv(mix, x)
%LOGGMMACTIV Computes the log of activations of a Gaussian mixture model.
%
%	Description
%	This function computes the activations A (i.e. the  probability
%	P(X|J) of the data conditioned on each component density)  for a
%	Gaussian mixture model.  For the PPCA model, each activation is the
%	conditional probability of X given that it is generated by the
%	component subspace. The data structure MIX defines the mixture model,
%	while the matrix X contains the data vectors.  Each row of X
%	represents a single vector.
%
%	See also
%	GMM, LOGGMMPOST, LOGGMMPROB
%

%	Copyright (c) Ian T Nabney (1996-9)
%       log implementation Copyright (c) Alexander G Dimitrov 2001

% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

ndata = size(x, 1);
a = zeros(ndata, mix.ncentres);  % Preallocate matrix

switch mix.covar_type
  
case 'spherical'
  % Calculate squared norm matrix, of dimension (ndata, ncentres)
  n2 = dist2(x, mix.centres);
  
  % Calculate width factors
  wi2 = ones(ndata, 1) * (2 .* mix.covars);
  normal = (pi .* wi2) .^ (mix.nin/2);
  
  % Now compute the activations
  a = (-(n2./wi2))- log(normal);
  
case 'diag'
  normal = (2*pi)^(mix.nin/2);
  s = prod(sqrt(mix.covars), 2);
  for i = 1:mix.ncentres
    diffs = x - (ones(ndata, 1) * mix.centres(i, :));
    a(:, i) = (-0.5*sum((diffs.*diffs)./(ones(ndata, 1) * ...
      mix.covars(i,:)), 2)) -log(normal*s(i));
  end
  
case 'full'
  log_normal = log(2*pi)*(mix.nin/2);
  for i = 1:mix.ncentres
    diffs = x - (ones(ndata, 1) * mix.centres(i,:));
    % Use Cholesky decomposition of covariance matrix to speed computation
    c = chol(mix.covars(:,:,i));
    temp = diffs/c;
    a(:,i) = sum(temp.*temp, 2)/(-2)-log_normal-sum(log(diag(c)));
  end
case 'ppca'
  log_normal = mix.nin*log(2*pi);
  d2 = zeros(ndata, mix.ncentres);
  logZ = zeros(1, mix.ncentres);
  for i = 1:mix.ncentres
    k = 1 - mix.covars(i)./mix.lambda(i, :);
    logZ(i) = log_normal + mix.nin*log(mix.covars(i)) - ...
      sum(log(1 - k));
    diffs = x - ones(ndata, 1)*mix.centres(i, :);
    proj = diffs*mix.U(:, :, i);
    d2(:,i) = (sum(diffs.*diffs, 2) - ...
      sum((proj.*(ones(ndata, 1)*k)).*proj, 2)) / ...
      mix.covars(i);
  end
  a = (d2 + ones(ndata, 1)*logZ)/(-2);
otherwise
  error(['Unknown covariance type ', mix.covar_type]);
end
  
