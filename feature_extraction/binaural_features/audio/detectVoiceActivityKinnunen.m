function vad = detectVoiceActivityKinnunen(in,thresdB,preMethod,format,blockSize,hopSize)
%detectVoiceActivityKinnunen   Energy-based voice activity detection.
%
%USAGE 
%         vad = detectVoiceActivityKinnunen(in)
%         vad = detectVoiceActivityKinnunen(in,thresdB,preMethod,format)
%
%INPUT ARGUMENTS
%          in : audio object (see readAudio or genAudio)
%     thresdB : energy threshold in dB
%   preMethod : string/cell defining pre-processing (see preprocessAudio)
%      format : output format of VAD decision ('samples' or 'frames')
%
%OUTPUT ARGUMENTS
%         vad : voice activity decision [nSamples|nFrames x 1]
%
%NOTE
%   If no output is specified, the VAD detection will be plotted.
% 
%REFERENCES
%   [1] T. Kinnunenand H. Lib, "An Overview of Text-Independent Speaker
%       Recognition: from Features to Supervectors",Speech Communication,
%       Vol.52, Issue 1, pp.12-40, 2010.

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2009 
%              TUe Eindhoven  
%              t.may@tue.nl   
%
%   History :
%   v.0.1   2009/10/31
%   v.0.2   2009/11/07 added VAD output format to input parameters
%   v.0.3   2012/05/23 prevent scaling of VAD decision in OLA framework
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
narginchk(1, 6);

% Set default values
if nargin < 6 || isempty(hopSize);    hopSize   = 10e-3;      end
if nargin < 5 || isempty(blockSize);  blockSize = 20e-3;      end
if nargin < 4 || isempty(format);     format    = 'frames';   end
if nargin < 3 || isempty(preMethod);  preMethod = 'removedc'; end
if nargin < 2 || isempty(thresdB);    thresdB   = 40;         end


%% **************************  CHECK AUDIO DATA  **************************
% 
% 
% Check if IN is an audio structure
if ~isAudio(in);
   error('Input argument must be an audio structure.') 
end

% Mono downmix
if in.nChannels > 1
    in = mono(in);
end


%% ****************************  VAD PARAMETER  ***************************
% 
% 
% Block processing
winType = 'rectwin';
nfft    = 'blockSize';

% Noise floor
noiseFloor = -55;

% Initialize framing parameter
P = audio2framesInit(in.fs,'blockSize',blockSize,'hopSize',hopSize,...
                     'winType',winType,'nfft',nfft);

% Prevent scaling of VAD decision                 
P.WOLA.reScale = 1;                 


%% ************************  DETECT VOICE ACTIVITY  ***********************
% 
% 
% Preprocessing
in = preprocessAudio(in,preMethod);

% Frame audio
data = audio2frames(in,P);

% Estimate energy in dB
energy = calcEnergy(data + eps,'decibel');

% Set max to 0dB
energy = energy - max(energy);
        
% VAD decision (frame-based)
frameVAD = energy > -abs(thresdB) & energy > noiseFloor;


%% *************************  RETURN VAD DECISION  ************************
% 
% 
% Select output format
switch lower(format)
    case 'samples'
        % Convert frame-based VAD decision to samples
        vad = frames2audio(repmat(frameVAD, [P.WOLA.blockSize 1]),P);
        
        % Quantize VAD decision (remove effect of analysis window)
        vad = round(vad);
        
        % Extend VAD data to match input size
        vad = logical([vad;  repmat(vad(end),[in.nSamples-length(vad) 1])]);
    case 'frames'
        % Return logical VAD decision
        vad = logical(frameVAD).';
    otherwise
        error(['VAD format "',lower(format),'" is not recognized.'])
end


%% ****************************  SHOW RESULTS  ****************************
% 
% 
% Plot results
if nargout == 0
    t    = (1:in.nSamples)/in.fs;
    tIdx = genAxisTime(P.WOLA.blockSize,P.WOLA.hopSize,...
                       length(frameVAD),in.fs);
   
    if isequal(lower(format),'frames')
        vIdx = tIdx;
    else
        vIdx = t;
    end
                              
    figure;
    ax(1) = subplot(3,1,[1 2]);
    plot(t,in.data)
    hold on;plot(vIdx,max(abs(in.data))*vad,'k','LineWidth',2)
    ylabel('Amplitude')
    xlim([0 inf])

    ax(2) = subplot(3,1,3);
    h = plot(tIdx,energy,[tIdx(1) tIdx(end)],-abs([thresdB thresdB]));
    set(h(2),'Color','k','LineStyle','--')
    legend({'energy' 'threshold'},'Location','SouthEast',...
           'Orientation','horizontal','FontSize',8)
    xlabel('Time (sec)')
    ylabel('Amplitude')
    ylim([-60 0])
    
    linkaxes(ax,'x');
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