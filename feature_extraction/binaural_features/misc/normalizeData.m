function [out,scale] = normalizeData(in,method)
%normalizeData   Perform channel-dependent normalization of input data.
%
%USAGE
%   [OUT,SCALE] = normalizeData(IN)
%   [OUT,SCALE] = normalizeData(IN,METHOD)
%   
%INPUT ARGUMENTS
%       IN : data matrix arranged as [nPoints x nChannels]
%   METHOD : string specifying normalization method
%            'mean'    - normalize data to have zero mean
%            'var'     - normalize data to have unit variance
%            'meanvar' - normalize data to have zero mean and unit variance
%            'max'     - normalize data to its maximum
% 
%            (default, METHOD = 'meanvar')
% 
%OUTPUT ARGUMENTS
%      OUT : normalized data [nPoints x nChannels]
%    SCALE : normalization factors per channel [1|2 x nChannels]
% 
%NOTE
%   The normalization factor is computed across all "nPoints" for each
%   channel independently. 

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2009 
%              TUe Eindhoven  
%              t.may@tue.nl   
%
%   History :
%   v.0.1   2009/10/12
%   ***********************************************************************

% Check for proper input arguments
narginchk(1, 2);

% Set default value
if nargin < 2 ||isempty(method); method = 'meanvar'; end;

% Select normalization method
switch lower(method)
    case 'mean'
        % Zero mean
        [out,scale] = normalizeMEX(in,1);
    case 'var'
        % Unit variance
        [out,scale] = normalizeMEX(in,2);
    case 'meanvar'
        % Zero mean and unit variance
        [out,scale] = normalizeMEX(in,3);
    case 'max'
        % Maximum value of one
        scale = max(abs(in),[],1);
        
        % Check normalization constant
        if any(scale==0) || any(~isfinite(scale))
            error('Normalization constant is zero/not finite')
        end
        
        % Perform normalization
        out = in./repmat(scale,[size(in,1) 1]);
    otherwise
        error(['Normalization method "',lower(method),'" is not ',...
               'recognized']);
end 


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