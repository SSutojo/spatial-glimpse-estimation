function fObj = genFilterObj(b,a,varargin)
%genFilterObj   Initialize filter object.
% 
%USAGE
%   f = genFilterObj(b)
%   f = genFilterObj(b,a)
%   f = genFilterObj(b,a,PARNAME1,PARVAL1,PARNAME2,PARVAL2,...)
%
%INPUT ARGUMENTS
%   b : filter coefficients of the numerator
%   a : filter coefficients of the denominator (default, a = 1)
% 
%   The following PARNAME/PARVALUE pairs are recognized
%
%    PARNAME  | PARVALUE
%    ------------------------------------------
%        'fs' : sampling frequency in Hz
%     'label' : filter description (string)
% 
%OUTPUT ARGUMENTS
%   f : filter object
% 
%NOTE
%   The filter coefficients will be normalized to ensure that a(1) = 1.
%   Furthermore, it will be checked if the poles of the filter denominator
%   are within the unit circle to ensure stability. 
% 
%EXAMPLES
%   % Create first order difference filter
%   f = genFilterObj([1 -1],1,'label','first-order differentiator');
% 
%   % Plot filter response
%   plotFilter(f);
% 
%   See also plotFilter, calcResponseFilter, and cascadeFilter.


%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2007-2009 
%              TUe Eindhoven
%              t.may@tue.nl    
%
%   History :  
%   v.0.1   2007/08/05
%   v.0.2   2009/09/26
%   v.0.3   2009/10/15 cleaned up
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
narginchk(1, inf);

% Set default value
if nargin < 2 || isempty(a); a = 1; end

% Check if filter coefficients are vectors
if ~isvector(b) || ~isvector(a)
    error('Filter coefficients must be vectors.')
else
    % Ensure column vectors
    b = b(:); a = a(:);
end

% Check for proper denominator
if a(1) == 0
    error('First element of denominator must be non-zero.')
end

% Normalization
if a(1) ~= 1
    % Display warning
    warning('MATLAB:normalizeFilter',...
            'Normalize filter coefficients to enfore a(1) = 1.')
    
    % Rescale filter coefficients
    b = b/a(1);
    a = a/a(1);
end

% Determine filter type (FIR or IIR)
if length(a) > 1
    % Check for IIR filter stability
    if ~isStable(b,a)
        error('Filter poles are outside the unit cirlce.')
    end
    fType = 'IIR';
else
    fType = 'FIR';
end


%% **********************  CREATE FILTER STRUCTURE  ***********************

% Determine filter order
order = max(length(b),length(a))-1;

% Create filter states
states = zeros(order,1);                     

% Create filter object
fObj = struct('label','filter object','structure',...
              'Direct-Form II Transposed','type',fType,...
              'b',b,'a',a,'fs',false,...
              'order',order,'states',states,'cascade','');

% Parse additional input parameters
%
if ~isempty(varargin)
    % Prevent fields from beeing changed
    block = {'b' 'a' 'type' 'order' 'states' 'cascade'};
    
    % Parse input arguments
    fObj = parseInput2Struct(fObj,varargin,block);
end               


%% ****************************  PLOT FILTER  *****************************
if nargout == 0
    % Plot filter object
    plotFilter(fObj);
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