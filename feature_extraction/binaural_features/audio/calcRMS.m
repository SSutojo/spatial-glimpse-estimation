function out = calcRMS(input,dim) 
%calcRMS   Compute the Root-Mean-Square (RMS) value.
%
%USAGE 
%   OUT = calcRMS(INPUT)
%   OUT = calcRMS(INPUT,DIM)
%
%INPUT ARGUMENTS
%   INPUT : data matrix
%     DIM : DIM specifies the dimension along which the rms should be 
%           computed. If DIM is empty or omitted, the first non-singleton 
%           dimension is used for RMS calculation.
%
%OUTPUT PARAMETERS
%     OUT : RMS value of INPUT
% 
%   See also normalizeData.

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2009 
%              TUe Eindhoven  
%              t.may@tue.nl  
%
%   History :
%   v.0.1   2009/10/20
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
narginchk(1, 2);

% Set default value
if nargin < 2 || isempty(dim);  dim = find(size(input) ~= 1, 1); end

% Check dimensionality
if isempty(dim)
    error('INPUT is a scalar!');
end  
        
% Calculate RMS
out = sqrt(mean(input.^2,dim));


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