function bFilterObj = isFilterObj(varargin)
%isFilterObj   Check if input is a filter object. 
%
%USAGE 
%   bFilterObj = isFilterObj(in)
%
%INPUT ARGUMENTS
%   in : input parameter(s)
%
%OUTPUT ARGUMENTS
%   bFilterObj : true/false depending on whether in is a filter object.
%
%EXAMPLES
%   % Create filter object
%   fObj = genFilterRemoveDC(16e3);
%   % Check if "fObj" is a filter object
%   isFilterObj(fObj)
%   ans = 
%         1
%
%   isFilterObj(fObj,true)
%   ans = 
%         1
%         0
% 
%   See also genFilterObj.

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
narginchk(1, inf);


%% *****************  CHECK IF INPUT IS A FILTER OBJECT  ******************

% Required structure fields to be recognized as a filter object
reqFields = {'label','structure','type','b','a','fs',...
             'order','states','cascade'};

% One structure containing multiple filter objects
if length(varargin{1}) > 1 && isequal(length(varargin),1)
    if isstruct(varargin{1})
        % Number of input arguments
        nArgs = length(varargin{1});
        
        % Initialize output
        bFilterObj = false(nArgs,1);
        
        % Loop over number of input arguments
        for ii = 1 : nArgs
            % Check if IN is a structure
            if isstruct(varargin{1}(ii))
                
                % Check if all required fields are present
                if all(isfield(varargin{1}(ii),reqFields))
                    % Set flag to true
                    bFilterObj(ii) = true;
                end
            end
        end
    elseif iscell(varargin{1})
        % Number of input arguments
        nArgs = length(varargin{1});
        
        % Initialize output
        bFilterObj = false(nArgs,1);
        
        % Loop over number of input arguments
        for ii = 1 : nArgs
            % Check if IN is a structure
            if isstruct(varargin{1}{ii})
               
                % Check if all required fields are present
                if all(isfield(varargin{1}{ii},reqFields))
                    % Set flag to true
                    bFilterObj(ii) = true;
                end
            end
        end
    else
        error('Format of filter objects is not recognized.')
    end
else
    % Number of input arguments
    nArgs = length(varargin);

    % Initialize output
    bFilterObj = false(nArgs,1);

    % Loop over number of input arguments
    for ii = 1 : nArgs
        % Check if IN is a structure
        if isstruct(varargin{ii})
            
            % Check if all required fields are present
            if all(isfield(varargin{ii},reqFields))
                % Set flag to true
                bFilterObj(ii) = true;
            end
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