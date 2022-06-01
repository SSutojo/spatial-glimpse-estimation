function [mix,speech,noise,al] = controlSNR(speech,noise,snrdB,bVAD,weight)
%controlSNR   Mix audio signals at specified Signal-to-Noise Ratio (SNR).
%
%USAGE
%   [MIX,SPEECH,NOISE] = controlSNR(SPEECH,NOISE,SNRdB,bVAD,WEIGHT);
%
%INPUT ARGUMENTS
%      SPEECH : speech audio structure or filename
%       NOISE : noise audio structure or filename which should be added to
%               SPEECH. The level of NOISE is adapted according to the 
%               specified SNR, relative to SPEECH.
%       SNRdB : Signal-to-Noise Ratio (SNR) in decibels
%        bVAD : if true, reject silent signal parts for SNR calculation
%               if false, use complete signals for SNR calculation
%               The VAD can be set for both SPEECH and NOISE independently.
%               Per default, the VAD is only used for the SPEECH signal and
%               not for the noise. (default, bVAD = [true false])
%      WEIGHT : spectral weighting for SNR calculation
%               'a-weighted' 'b-weighted' 'c-weighted' 'd-weighted'
%               (default, WEIGHT = 'a-weighted')
% 
%OUTPUT ARGUMENTS
%         MIX : audio structure containing mixture of SPEECH and NOISE
%      SPEECH : speech audio structure
%       NOISE : noise audio structure with adapted level
%
%NOTE
%   The signal level of NOISE is changed according to the specified SNR,
%   whereas the level of SPEECH remains unchanged.


%   Developed with Matlab 7.5.0.342 (R2007b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2007-2010 
%              University of Oldenburg and TU/e Eindhoven 
%              tobias.may@uni-oldenburg.de   t.may@tue.nl
%
%   History :
%   v.0.1   2007/11/01
%   v.0.2   2008/03/13 suppress warning if file length differs
%   v.0.3   2008/05/17 clean-up
%   v.0.4   2008/07/01 added multi-channel SNR calculation
%   v.0.5   2008/12/11 added flag to truncate mixture if filelength differs
%   v.0.6   2009/05/11 added support for -inf and inf SNR
%   v.1.0   2009/08/07 added VAD support
%   v.1.1   2009/11/07 added SNR weighting
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
narginchk(3 ,5);

% Set default values
if nargin < 5 || isempty(weight); weight = 'a-weighted'; end
if nargin < 4 || isempty(bVAD);   bVAD   = [true false]; end


%% ***********************  CHECK AUDIO SIGNALS  **************************
% 
% 
% Check if first two input variables are strings/cells
if (ischar(speech) && ischar(noise)) || (iscell(speech) && iscell(noise))
    % Read audio files
    speech = readAudio(speech);
    noise  = readAudio(noise);
elseif any([~isAudio(speech) ~isAudio(noise)])
    % Check if input is a struct 
    error(['The first two input arguments must be either two strings',...
           ' specifying wave files or two audio structures.'])
else
    % Update audio parameters
    speech = updateAudio(speech);
    noise  = updateAudio(noise);
end

% Check for coherent signal characteristics
if speech.fs ~= noise.fs
    error('Sampling rate mismatch.');
end

% Check number of bits
if isfield(speech,'nBits') && ...
   isfield(noise,'nBits')  && ...
   speech.nBits ~= noise.nBits
       error('Number of bits mismatch!');
end

% Flag for reswapping ...
bReSwap = false;

% Check SNR parameter
if isequal(snrdB,inf)
    % Return speech only
    bCompSNR = false;
elseif isequal(snrdB,-inf)
    % Return noise only, therefore swap speech and noise
    bCompSNR = false;
    % Set flat
    bReSwap = true;
    % Swap speech and noise
    [speech,noise] = swap(speech,noise);
elseif ~isnumeric(snrdB)
    error('SNR must be numeric.');
else
    % Compute SNR
    bCompSNR = true;
end


%% *********************  APPLY SPECTRAL WEIGHTING  ***********************
%
% 
% Filter speech and noise signals
speechW = preprocessAudio(speech,weight);
noiseW  = preprocessAudio(noise,weight);


%% ***************************  COMPUTE SNR  ******************************
% 
% 
% Initialize VAD 
switch length(bVAD)
    case 1
        bVAD = repmat(bVAD,[1 2]);
    case 2
        % Decide usage of VAD independently for speech and noise
    otherwise
        error(['VAD parameter must contain either 1 or 2 ',...
               'logical elements defining whether a VAD ',...
               'should be applied prior to SNR calculation.'])
end

% Compute SNR
if bCompSNR
    % Initialize VAD parameter
    if any(bVAD)
        % VAD threshold
        thresdB = 40;
        % Pre-process
        preMethod = ''; 
        % VAD output format
        f = 'samples';
    end
    
    % Perform Voice Activity Detection
    if bVAD(1)
        % VAD for speech
        vadS = detectVoiceActivityKinnunen(speech,thresdB,preMethod,f);
    else
        % Use all samples
        vadS = true(speech.nSamples,1);
    end
    if bVAD(2)
        % VAD for noise
        vadN = detectVoiceActivityKinnunen(noise,thresdB,preMethod,f);
    else
        % Use all samples
        vadN = true(noise.nSamples,1);
    end

    % Get average SNR over all channels by summing the calculated energy
    % over all channels => Heidi Christensen, Interspeech paper 2010
    currentSNR = sqrt(sum(mean(power(speechW.data(vadS,:),2),1),2) ./ ...
                      sum(mean(power(noiseW.data(vadN,:),2),1),2));
             
    % Transform SNR in decibel to linear value
    snrFactor = 10^(-snrdB/20);
else
    % Used if SNR = inf, SNR = -inf  
    currentSNR = 0;
    snrFactor  = 0;
end


%% ****************************  ADJUST SNR  ******************************
% 
% 
al = currentSNR * snrFactor;

% Change noise level according to current SNR and the specified SNR
noise.data = noise.data * al;

% Add signals
mix = addAudio(speech,noise);

% Update audio label
mix.label = ['mix@',num2str(snrdB),'dB'];

% Re-swap speech and noise
if bReSwap
    [speech,noise] = swap(speech,noise);
end


%% ***************************  PLOT MIXTURE  *****************************
% 
% 
% Plot mixture
if nargout == 0
   plotAudio(mix);
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