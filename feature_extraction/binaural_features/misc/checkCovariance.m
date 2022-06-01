function ok = checkCovariance(M)



%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
narginchk(1, 1);

% Check for proper input 
if ~all(isfield(M,{'ncentres' 'covars'}))
   error('Input must be a GMM NETLAB structure') 
end


%% **********************  CHECK COVARIANCE MATRIX  ***********************

% Allocate memory
ok = zeros(M.ncentres,1);

% Checks if the covariance of a NETLAB GMM structure M are ok (i.e., do not
% result in a singular or close to singular covariance matrix)
for ii = 1 : M.ncentres
    covs = M.covars(ii,:);
    C    = diag(covs);
    % LAPACK reciprocal condition estimator
    if ~isnan(rcond(C))
      ok(ii) = 1;
    else
        % Set weight to zero
        M.priors(ii) = 0;
        M.priors = M.priors ./ sum(M.priors);
    end
end

% Accumulate results over all Gaussian components
ok = prod(ok);


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************
