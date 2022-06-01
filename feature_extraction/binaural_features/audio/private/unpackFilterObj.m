function fObj = unpackFilterObj(fObj)
%unpackFilterObj   Unpack and check filter objects.
%
%USAGE 
%   fObj = unpackFilterObj(fObj)
%
%INPUT ARGUMENTS
%   fObj : structure/cell containing filter objects
%
%OUTPUT ARGUMENTS
%   fObj : structure containing filter objects
%
%   See also ...

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
narginchk(1, 1);


%% ****************  UNPACK AND CHECK ALL FILTER OBJECTS  *****************

% Unpack all filter objects
if isstruct(fObj)
    % ... proper format
    
    % Check if all "fObj" elements are filter objects
    if ~all(isFilterObj(fObj))
        error('Not all input arguments are filter objects.')
    end
    
elseif iscell(fObj)
    % Check if all "fObj" elements are filter objects
    if ~all(isFilterObj(fObj{:}))
        error('Not all input arguments are filter objects.')
    else
        % Unpack filter objects
        if isstruct(fObj{1})
            fObj = [fObj{:}];
        elseif iscell(fObj{1})
            fObj = [fObj{:}{:}];
        else
            error(['Format of filter objects is not recognized.'])
        end
    end    
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