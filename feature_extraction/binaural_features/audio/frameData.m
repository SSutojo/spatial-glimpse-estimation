function out = frameData(input,blockSize,hopSize,win)
%frameData   Frame data.
% 
%USAGE
%   frames = frameData(input,blockSize,hopSize,window,nFrames);
% 
%INPUT ARGUMENTS
% 
%OUTPUT ARGUMENTS
% 
%EXAMPLE
% 

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2009 
%              TUe Eindhoven  
%              t.may@tue.nl     
%
%   History :
%   v.0.1   2009/11/23
%   v.0.2   2010/05/18 automatically determine the number of frames
%   ***********************************************************************

% Check for proper input arguments
narginchk(4, 4);

% Compute number of frames
nFrames = max(fix((size(input,1)-(blockSize-hopSize))/(hopSize)),1);

% Check if window is a window in samples or a string
if ischar(win)
    win = window(win,blockSize);
else
    if blockSize ~= length(win)    
       error('Mismatch between blocksize and window size.') 
    end
end

% Framing (MEX processing)
out = frameDataMEX(input,blockSize,hopSize,win,nFrames);


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