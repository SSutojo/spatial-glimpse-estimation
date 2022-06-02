function [post, a] = loggmmpost(mix, x)
%function [post, a] = loggmmpost(mix, x)
% logGMMPOST Computes the log of class posterior probabilities of a
% Gaussian mixture model. 
%
%	Description
%	This function computes the posteriors POST (i.e. the probability of
%	each component conditioned on the data P(J|X)) for a Gaussian mixture
%	model.   The data structure MIX defines the mixture model, while the
%	matrix X contains the data vectors.  Each row of X represents a
%	single vector.
%
%	See also
%	logGMM, logGMMACTIV, logGMMPROB
%

%	Copyright (c) Ian T Nabney (1996-9)
%       log implementation Copyright (c) Alexander G Dimitrov 2001

% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

ndata = size(x, 1);

a = loggmmactiv(mix, x);

post = log(ones(ndata, 1)*mix.priors) + a;
Maxpost=max(post,[],2); post=post-repmat(Maxpost,1,size(post,2));
% renormalize here, make max post ~ 1
s = sum(exp(post), 2); 
% AGD: Set any zeros to one before dividing. If there is a zero in s,
% then the whole row of post should be zero. Not very likely here,
% because of the logs, but keep it just in case. Likely to misfire at
% times.
s = s + (s==0);
s=log(s);
% I don't have to correct for Maxpost here, since in the last line,
% the two Maxposts for post and s cancel out.
post = post-s*ones(1, mix.ncentres);
