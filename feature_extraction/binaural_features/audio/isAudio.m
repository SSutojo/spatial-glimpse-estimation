function out = isAudio(in)
%isAudio   Check if input is an audio object. 
%   An audio object is a MATLAB structure with at least two fiels, namely
%   "DATA" and "FS". 
%
%USAGE
%   OUT = isAudio(IN)
%   
%INPUT ARGUMENTS
%    IN : input
% 
%OUTPUT ARGUMENTS
%   OUT : true/false depending on whether IN is an audio structure.
% 
%EXAMPLE
%   % Load audio 
%   load handel
%   % Check 
%   isAudio(struct('data',y,'fs',Fs))
%   ans = 
%         1
% 
%   See also genAudio, readAudio and plotAudio.

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2008-2009
%              TUe Eindhoven 
%              t.may@tue.nl    
% 
%   History :  
%   v.0.1   2008/05/14
%   v.0.2   2009/10/20
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
narginchk(1, 1);


%% *****************  CHECK IF IN IS AN AUDIO STRUCTURE  ******************

% Initialize output
out = false; 

% Check if IN is a structure
if isstruct(in)
    % Required structure fields
    reqFields = {'data' 'fs'};
    
    % Check if all required fields are present
    if all(isfield(in,reqFields))
        % Set flag to true
        out = true;
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