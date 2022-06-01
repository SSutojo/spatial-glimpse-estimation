function timeSec = genAxisTime(blockSize,hopSize,nBlocks,fs)
%genAxisTime   Create time axis in seconds.
%
%USAGE
%     timeSec = genAxisTime(blockSize,hopSize,nBlocks,fs)
%
%INPUT ARGUMENTS
%   blockSize : block size in samples
%     hopSize : hop size in samples
%     nBlocks : number of frames
%          fs : sampling frequency in Hz
% 
%OUTPUT ARGUMENTS
%     timeSec : time axis in seconds
% 
%EXAMPLE
%   timeSec = genAxisTime(320,160,150,16e3);
% 
%   See also genAxisFreq.

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2008-2011 
%              University of Oldenburg and TU/e Eindhoven   
%              tobias.may@uni-oldenburg.de   t.may@tue.nl
%
%   History :  
%   v.0.1   2008/05/11
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
narginchk(4,4);


%% *************************  CREATE TIME AXIS  ***************************
% 
% 
% Block axis
colIdx = (0:(nBlocks-1))*(hopSize);

% Time vector whose elements are centered 
timeSec = (colIdx + ((blockSize)/2))/fs; 


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
