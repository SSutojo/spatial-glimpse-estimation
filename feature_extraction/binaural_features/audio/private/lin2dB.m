function dB = lin2dB(in,method,dataFloor,ref)
%lin2dB   Convert power/magnitude values to decibel units.
%
%USAGE
%   dB = lin2dB(inLin)
%   dB = lin2dB(inLin,method,dataFloor,ref)
%
%INPUT ARGUMENTS
%       inLin : data matrix 
%      method : select conversion method
%               'magnitude' - magnitude to decibel conversion (default)
%                   'power' - power to decibel conversion 
%   dataFloor : floor data prior to dB conversion (default, floor = 0)
%         ref : reference (default, ref = 1)
% 
%OUTPUT ARGUMENTS
%       lin : data matrix in decibel units
% 
%   See also dB2lin.


%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2008-2009
%              TUe Eindhoven   
%              t.may@tue.nl   
%
%   History :  
%   v.0.1   2008/05/16
%   v.0.2   2009/10/21 cleaned up
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
narginchk(1, 4);

% Set default values
if nargin < 4 || isempty(ref);       ref       = 1;           end
if nargin < 3 || isempty(dataFloor); dataFloor = 0;           end
if nargin < 2 || isempty(method);    method    = 'magnitude'; end


%% ************************  CONVERT INPUT TO dB  *************************

% Set decibel floor
if dataFloor > 0 
    in = max(abs(in),dataFloor);
end

% Select method
switch lower(method)
    case 'magnitude'
        dB = 20 * log10(abs(in)./ref);
    case 'power'
        dB = 10 * log10(abs(in)./ref);
    otherwise
        error(['Decibel conversion "',method,'" is not recognized.'])
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