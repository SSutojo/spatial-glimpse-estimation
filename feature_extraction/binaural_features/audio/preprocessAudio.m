function audio = preprocessAudio(audio,preMethod)
%preprocessAudio   Preprocess audio object.
%
%USAGE 
%   audio = preprocessAudio(audio,preMethod)
%
%INPUT ARGUMENTS
%       audio : audio object
%   preMethod : string defining preprocessing method
%               '' or 'nothing' = do nothing
%               'normrms'       = normalize audio with respect to its RMS
%               'normmax'       = normalize audio with respect to its
%                                 maximum 
%               'whiteningfir'  = Apply first order whitening filter
%               'whiteningfft'  = Perform FFT-based whitening
%               'a-weighted'    = Apply A-weighting filter
%               'b-weighted'    = Apply B-weighting filter
%               'c-weighted'    = Apply C-weighting filter
%               'd-weighted'    = Apply D-weighting filter
%               'removedc'      = Use 4th order high-pass with a cutoff
%                                 frequency of 20 Hz to remove DC omponents
%               'lowpass'       = Apply 5th order low-pass with a cutoff
%                                 frequency of 4500 Hz
% 
%OUTPUT ARGUMENTS
%       audio : processed audio object
% 
%NOTE
%   It's possible to combine various preprocessing operations by using a
%   cell array of methods {'whiteningfft' 'a-weighted'}


%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2009 
%              TUe Eindhoven 
%              t.may@tue.nl     
%
%   History :
%   v.0.1   2009/09/025
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
narginchk(2,2);

% Check audio
if ischar(audio) || (iscell(audio))
    % Read audio file
    audio = readAudio(audio);
elseif ~isAudio(audio) 
    % Report error
    error(['''Audio'' must be either a string/cell',...
           ' specifying a wave file or an audio structure.'])
else
    % Update audio parameters
    audio = updateAudio(audio);
end

% Check if preMethod is char or cell array
if ischar(preMethod) || islogical(preMethod)
    preMethod = {preMethod};
    nMethods  = 1;
else
    nMethods  = length(preMethod);
end


%% *************************  PERFORM FILTERING  **************************
% 
% 
% Loop over number of preprocessing methods
for ii = 1 : nMethods
    % Select preprocessing method
    switch lower(preMethod{ii})
        case {'' 'nothing' false}
            % Do nothing ...
            
        case 'normrms'
            % Normalize audio signal with respect to its rms 
            audio = normalizeAudio(audio,'rms');
        case 'normmax'
            % Normalize audio signal with respect to its maximum
            audio = normalizeAudio(audio,'max');
        case 'whiteningfir'
            % Apply 1st order whitening filter
             audio.data = filter([1 -.97], 1, audio.data);
        case 'whiteningfft'
            % Channel-dependent whitening
            for jj = 1 : audio.nChannels
                spec             = fft(audio.data(:,jj));
                audio.data(:,jj) = real(ifft(spec./abs(spec)));
            end
        case {'a-weighted' 'b-weighted' 'c-weighted' 'd-weighted'}
            % Create A,B,C or D filter
            F = genFilterABCD(audio.fs,preMethod{ii}(1));

            % Apply weighting filter
            audio = filterAudio(audio,F);
        case 'iem'
            % Generate middle ear filter
            F = genFilterMiddleEar(audio.fs,'IEM');
            
            % Filtering
            audio = filterAudio(audio,F);
        case 'removedc'
            % Cutoff frequency in Hertz
            cf = 20;
            
            % Create high-pass filter
            F = genFilterRemoveDC(audio.fs,cf);
            
            % Apply highpass filter
            audio = filterAudio(audio,F);

        case 'lowpass'
            % Cutoff frequency in Hertz
            cutoffHz = 5000;
            
            % Design 5th order butterworth lowpass
            [b, a] = butter(5, cutoffHz / (0.5 * audio.fs),'low');
            
            % Create filter object
            F = genFilterObj(b,a,'fs',audio.fs,'label','low-pass filter');
            
            % Apply highpass filter
            audio = filterAudio(audio,F);
            
        otherwise
            error(['Pre-processing method ''',lower(preMethod{ii}),...
                   ''' does not exist.'])
    end
end


%% *************************  PLOT AUDIO SIGNAL  **************************
% 
% 
if nargout == 0
   plotAudio(audio);
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