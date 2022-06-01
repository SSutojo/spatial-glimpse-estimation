function [mix, options, errlog] = loggmmem(mix, x, options)
%LOGGMMEM	EM algorithm for Gaussian mixture model.
%
%	Description
%	[MIX, OPTIONS, ERRLOG] = GMMEM(MIX, X, OPTIONS) uses the Expectation
%	Maximization algorithm of Dempster et al. to estimate the parameters
%	of a Gaussian mixture model defined by a data structure MIX. The
%	matrix X represents the data whose expectation is maximized, with
%	each row corresponding to a vector.    The optional parameters have
%	the following interpretations.
%
%	OPTIONS(1) is set to 1 to display error values; also logs error
%	values in the return argument ERRLOG. If OPTIONS(1) is set to 0, then
%	only warning messages are displayed.  If OPTIONS(1) is -1, then
%	nothing is displayed.
%
%	OPTIONS(3) is a measure of the absolute precision required of the
%	error function at the solution. If the change in log likelihood
%	between two steps of the EM algorithm is less than this value, then
%	the function terminates.
%
%	OPTIONS(5) is set to 1 if a covariance matrix is reset to its
%	original value when any of its singular values are too small (less
%	than MIN_COVAR which has the value eps).   With the default value of
%	0 no action is taken.
%
%	OPTIONS(14) is the maximum number of iterations; default 100.
%
%	The optional return value OPTIONS contains the final error value
%	(i.e. data log likelihood) in OPTIONS(8).
%
%	See also
%	GMM, GMMINIT
%

%	Copyright (c) Ian T Nabney (1996-9)
%       log implementation Copyright (c) Alexander G Dimitrov 2001

% log version is useful for high-dimensional models, where many of the
% parameters in the original gmmem hit numeric zeros.  tested against
% gmmem with a 4 element 2d full model (logtest.mat). Practically
% indistinguishable (abs (em_err - logem_err)<1e-15); length(em_err) ~
% 100.
% 
% error here is given per sample, so different dataset sizes
% affect the precision less. maybe rescaling per dimension needed
% as well.

% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

[ndata, xdim] = size(x);

% Sort out the options
if (options(14))
  niters = options(14);
else
  niters = 100;
end

display = options(1);
store = 0;
if (nargout > 2)
  store = 1;	% Store the error values to return them
  errlog = zeros(1, niters);
end
test = 0;
if options(3) > 0.0
  test = 1;	% Test log likelihood for termination
end

check_covars = 0;
if options(5) >= 1
  if display >= 0
    disp('check_covars is on');
  end
  check_covars = 1;	% Ensure that covariances don't collapse
  MIN_COVAR = eps;	% Minimum singular value of covariance matrix
  init_covars = mix.covars;
end

% Main loop of algorithm. Differences from gmmem start here.
for n = 1:niters
  
  % Calculate posteriors based on old parameters
%  [post, act] = gmmpost(mix, x);
   [logpost, act] = loggmmpost(mix, x);

  % Calculate error value if needed
  if (display | store | test)
    logprob = logtimes(act,mix.priors');
    % Error value is negative log likelihood of data
    e = - sum(logprob)/ndata;  % per sample
    if store
      errlog(n) = e;
    end
    if display > 0
      fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
    end
    if test
      if (n > 1 & abs(e - eold) < options(3))
        options(8) = e;
        return;
      else
        eold = e;
      end
    end
  end
  
  % Adjust the new estimates for the parameters
  Maxpost=max(logpost(:)); %logpost=logpost-Maxpost;
  %new_pr = sum(exp(logpost), 1).*exp(Maxpost);
  lognew_pr=  logtimes(logpost',ones(size(logpost,1),1));  %log(sum(post,1))
  lognew_c =  logtimes(logpost', x);
  
  % Now move new estimates to old parameter vectors
  mix.priors  = rexp( lognew_pr' - log(ndata) );  % keep priors, not log priors
  mix.centres = rexp( lognew_c - repmat(lognew_pr,[1 mix.nin]) );
  
  switch mix.covar_type
  case 'spherical'
    n2 = dist2(x, mix.centres);
    for j = 1:mix.ncentres
      v(j) = logtimes(logpost(:,j)',n2(:,j)); % use logpost.
    end
    mix.covars = rexp((v-lognew_pr')-log(mix.nin));
    if check_covars
      % Ensure that no covariance is too small
      for j = 1:mix.ncentres
        if mix.covars(j) < MIN_COVAR
          mix.covars(j) = init_covars(j);
        end
      end
    end
  case 'diag'
    for j = 1:mix.ncentres
      diffs = x - (ones(ndata, 1) * mix.centres(j,:));
      mix.covars(j,:) = sum((diffs.*diffs).*...   % use logpost
	    repmat(rexp(logpost(:,j)-lognew_pr(j)),1,mix.nin) );
    end
    if check_covars
      % Ensure that no covariance is too small
      for j = 1:mix.ncentres
        if min(mix.covars(j,:)) < MIN_COVAR
          mix.covars(j,:) = init_covars(j,:);
        end
      end
    end
  case 'full'
    for j = 1:mix.ncentres
      diffs = x - (ones(ndata, 1) * mix.centres(j,:));
      diffs = diffs.*repmat(rexp((logpost(:,j)-lognew_pr(j)-Maxpost)/2),...
			    1, mix.nin); %use logpost
      mix.covars(:,:,j) = (diffs'*diffs)*rexp(Maxpost);
      % RC
    end
    if check_covars
      % Ensure that no covariance is too small
      for j = 1:mix.ncentres
        if min(svd(mix.covars(:,:,j))) < MIN_COVAR
          mix.covars(:,:,j) = init_covars(:,:,j);
        end
      end
    end
  case 'ppca'
    for j = 1:mix.ncentres
      diffs = x - (ones(ndata, 1) * mix.centres(j,:));
      diffs = diffs.*repmat(rexp( (logpost(:,j)-lognew_pr(j)-Maxpost)/2),...
				 1, mix.nin); %use logpost
      [mix.covars(j), mix.U(:,:,j), mix.lambda(j,:)] = ...
        ppca( (diffs'*diffs)*rexp(Maxpost), mix.ppca_dim);
      % AGD: see how ppca works. maybe I can rescale by exp(Maxpost)
      % after ppca is done...
    end
    if check_covars
      if mix.covars(j) < MIN_COVAR
        mix.covars(j) = init_covars(j);
      end
    end
    otherwise
      error(['Unknown covariance type ', mix.covar_type]);               
  end
end

if nargout>1
  %options(8) = -sum(log(gmmprob(mix, x)));
  options(8) = -sum(loggmmprob(mix, x));
end

if (display >= 0)
  disp('Warning: Maximum number of iterations has been exceeded');
end


function b=rexp(a)
% function b=rexp(a) = real(exp(a))

b=real(exp(a));