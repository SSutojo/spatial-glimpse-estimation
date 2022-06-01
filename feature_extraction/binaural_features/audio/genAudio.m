function aObj = genAudio(data,fs,label)
%genAudio   Create an audio object.
%   An audio object is a MATLAB structure with at least two fiels, namely 
%   'DATA' and 'FS'. 
%   
%USAGE
%    aObj = genAudio(DATA,FS)
%    aObj = genAudio(DATA,FS,LABEL)
%
%INPUT ARGUMENTS
%    DATA : audio signal arranged as [nSampels x nChannels]
%      FS : sampling frequency in Hz
%   LABEL : label specifying the audio signal (default, LABEL = '')
%
%OUTPUT ARGUMENTS
%    aObj : audio structure containing the following fields:
%        .label     - label
%        .data      - audio signal [nSampels x nChannels]
%        .fs        - sampling frequency in Hz
%        .nBits     - number of bits
%        .nSamples  - number of samples
%        .nChannels - number of channels
%        .range     - sample range of audio signal 
%        .format    - format description
%
%NOTE
%   The number of bits is determined by reading out the class information
%   of "data". The format of the audio structure is always set to 'MATLAB
%   workspace'. 
%
%EXAMPLE
%   % Load audio material
%   load chirp
% 
%   % Create audio object
%   in = genAudio(y,Fs,'chirp')
% 
%   See also isAudio, updateAudio, readAudio, normalizeAudio, plotAudio 
%            and resampleAudio.

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May
%              TUe Eindhoven 
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :  
%   v.0.1   2008/05/16
%   v.0.2   2009/10/20 cleaned up
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
narginchk(2, 3);

% Set default value
if nargin < 3 || isempty(label); label = []; end


%% ************************  CREATE AUDIO OBJECT  *************************

% Create and update audio object 
aObj = updateAudio(struct('data',data,'fs',fs,...
                          'format','MATLAB workspace'),label);

              
%% **************************  PLOT AUDIO OBJECT  *************************

% Plot signal
if nargout == 0
    plotAudio(aObj);
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