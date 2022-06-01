function in = normalizeAudio(in,method,bSeparately)
%normalizeAudio   Normalize audio object.
%
%USAGE
%   aObj = normalizeAudio(aObj)
%   aObj = normalizeAudio(aObj,method,bSeparately)
%   
%INPUT ARGUMENTS
%          aObj : audio object 
%        method : specify normalization measure (default, method = 'rms')
%                 'max' - normalize audio with respect to its maximum, such
%                         that all values are in the range [-1,1].
%                 'rms' - normalize audio with respect to its RMS
%   bSeparately : logical flag specifying multi-channel processing
%                  true - normalize channels independently
%                 false - preserve inter-channel relation by normalizing
%                         all channels to the maximum measure (e.g. rms)
%                         observed across channels
%                         (default, bSeparately = false)
% 
%OUTPUT ARGUMENTS
%      aObj : normalized audio object
% 
%EXAMPLE
%   % Load audio 
%   load handel
%   % Create audio structure
%   in = genAudio(y,Fs);
%   % Normalize audio object with respect to its rms
%   out = normalizeAudio(in,'rms');
%   % Check 
%   calcRMS(out.data)
%   ans = 
%         1
% 
%   See also genAudio, isAudio, updateAudio and plotAudio.


%   Developed with Matlab 7.6.0.324 (R2008a). Please send bug reports to:
%   
%   Author  :  Tobias May, 2008-2009 
%              TUe Eindhoven 
%              t.may@tue.nl    
%
%   History :  
%   v.0.1   2008/11/17
%   v.0.2   2009/02/26 added RMS normalization
%   v.0.3   2009/10/20 cleaned up
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
narginchk(1, 3);

% Set default values
if nargin < 2 || isempty(method);      method      = 'rms'; end
if nargin < 3 || isempty(bSeparately); bSeparately = false; end

% Check if IN is an audio structure
if ~isAudio(in);
   error('Input argument must be an audio structure.') 
end


%% ******************  EXTRACT NORMALIZATION MEASURE  *********************

% Select normalization method
switch lower(method)
    case 'max'
        % Compute signal maximum 
        normFactor = max(abs(in.data),[],1);
    case 'rms'
        % Compute RMS
        normFactor = calcRMS(in.data,1);
    otherwise
        error(['Normalization method "',lower(method),...
               '" is not recognized.'])
end


%% ***********************  NORMALIZE AUDIO DATA  *************************
if bSeparately
    % Handle each channel separately
else
    % Handle all channels together
    normFactor = repmat(max(normFactor),[in.nChannels 1]);
end

% Check normalization constant
if any(normFactor==0) || any(~isfinite(normFactor))
   error('Normalization constant is zero/not finite') 
end

% Loop over number of channels
for ii = 1 : in.nChannels
    in.data(:,ii) = in.data(:,ii)/normFactor(ii);
end


%% **************************  PLOT AUDIO DATA  ***************************

% Plot signal
if nargout == 0
    plotAudio(in);
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