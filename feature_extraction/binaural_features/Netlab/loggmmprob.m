function logprob = loggmmprob(mix, x)
%logGMMPROB Computes the log of data probability for a Gaussian mixture model.
%
%	Description
%	 This function computes the unconditional data density P(X) for a
%	Gaussian mixture model.  The data structure MIX defines the mixture
%	model, while the matrix X contains the data vectors.  Each row of X
%	represents a single vector.
%
%	See also
%	logGMM, logGMMPOST, logGMMACTIV
%

%	Copyright (c) Ian T Nabney (1996-9)
%       log implementation Copyright (c) Alexander G Dimitrov 2001

% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

% Compute activations
a = loggmmactiv(mix, x);
% Form dot product with priors
logprob = logtimes(a,mix.priors');   % log \sum_j p(x|j) p(j)