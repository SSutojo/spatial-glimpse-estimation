function P = parseInput2Struct(P,inArg,blockFields)
%parseInput2Struct   Parse parameter/value pairs to a structure.
%   The "parameters" are used to find the names of the fields in "P" which
%   should be updated with the corresponding "values".
%
%USAGE
%   P = parseInput2Struct(P,INARG)
%   P = parseInput2Struct(P,INARG,blockFields)
%      
%INPUT ARGUMENTS
%             P : structure which fields should be configurated according
%                 to the parameter/value pairs. 
%         INARG : cell array consisting of parameter/value pairs 
%                 (e.g. obtained from MATLAB's varargin)
%   BLOCKFIELDS : list of parameter which should not be updated by the
%                 parameter/value pairs. 
%
%OUTPUT ARGUMENTS
%             P : updated structure
%
%EXAMPLES
%   % Create parameter structure
%   P = struct('method','erb','range',[0 4000],'nFilter',32);
% 
%   % Update parameter structure 
%   P = parseInput2Struct(P,{'range',[50 5000]})
%   P = parseInput2Struct(P,{'range',[50 5000],'method','mel'})

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2008-2009 
%              TUe Eindhoven 
%              t.may@tue.nl   
%
%   History :
%   v.0.1   2008/05/18
%   v.0.2   2009/10/15 cleaned up
%   ***********************************************************************

% Check for proper input arguments
narginchk(2, 3);

% Set default value
if nargin < 3 || isempty(blockFields); blockFields = ''; end

% Get number of input arguemtns to parse
nArg = length(inArg);

% Check for proper dimension
if rem(nArg,2) ~= 0
    error('Field and value input arguments must come in pairs.')
end

% Loop over number of input arguments
for ii = 1 : 2 : nArg
   % Check if field 'inArg{ii}' does exist ...
   if isfield(P,inArg{ii})    
       % Check if field name should not be changed ...
       if strcmp(blockFields,inArg{ii})
           error(['Field name "',inArg{ii},'" cannot be changed.'])
       else
           % Configure P
           P = setfield(P,inArg{ii},inArg{ii+1});
       end
   else
       % Field does not exist ... so give out a warning
       warning('MATLAB:parseInput2Struct',['Field name "',inArg{ii},...
               '" is not recognized.'])
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