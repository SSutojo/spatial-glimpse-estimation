function [BinauralRatio, BetterEarSNR, BinauralUnmaskingAdvantage, BESNR_perTFFB, BUAdv_perTFFB] = ...
    vicente2020(FS, fc, MaskerSig, TargetSig, InternalNoise_L, InternalNoise_R, ceiling, HannWindowBE, HannWindowBU, WindowOverlap, weightings)
%VICENTE2020 Compute the effective SNR taking into account BU and BE
%   Usage: [predicted_SNR, BE, BU] = vicente2020(target_in,int_in,fs)
%
%   Input parameters:
%     target_in     : target
%     int_in        : interferer
%     fs            : sampling frequency [Hz]
%
%   Output parameters:
%     predicted_SNR : SNR predicted by the model
%     BE            : better-ear advantage
%     BU            : binaural masking level difference advantage
%
%   VICENTE2020 computes the effective SNR taking into account audibility (audiogram),
%   BU and BE by respective time frames, taking the target and interferer signals 
%   (sampled at fs) as inputs, along with the internal noise parameters computed 
%   using vicente2020_internalnoise
%
%   See also: lavandier2022 vicente2020nh vicente2020 prudhomme2020 leclere2015
%   jelfs2011
%
%   References:
%     M. Lavandier. A series of speech intelligibility models in the auditory
%     modeling toolbox. actaunited, 2022.
%     
%     Collin and Lavandier. Binaural speech intelligibility in rooms with
%     variations in spatial location of sources and modulation depth of noise
%     interferers. jasa, 134(2):1146--1159, 2013.
%     
%     Vicente and Lavandier. Further validation of a binaural model
%     predicting speech intelligibility against envelope-modulated noises.
%     hearingresearch, 390(107937), 2020.
%     
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/models/vicente2020.php

% Copyright (C) 2009-2021 Piotr Majdak, Clara Hollomey, and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 1.1.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   #StatusDoc: Perfect
%   #StatusCode: Good
%   #Verification: Qualified
%   #Requirements: MATLAB

%   AUTHOR: Matthieu Lavandier
%   adapted for AMT by Clara Hollomey (2021)


% Compute target info on the long term signal
[TargetLeftSpectrum, TargetRightSpectrum, TargetInterauralPhase] = local_gettargetstats2(TargetSig,fc,FS);

% Compute Better-Ear SNR
L = length(HannWindowBE);
OverlapSample = WindowOverlap*length(HannWindowBE);
BufferedLeft = buffer(MaskerSig(:,1), L, OverlapSample, 'nodelay') .* HannWindowBE;
BufferedRight = buffer(MaskerSig(:,2), L, OverlapSample, 'nodelay') .* HannWindowBE;

weighted_BE_SNR_frames = zeros(size(BufferedLeft,2),1);
BESNR_perTFFB = zeros(size(BufferedLeft,2),length(fc));

for n = 1:size(BufferedLeft,2)
    [weighted_BE_SNR_frames(n),BESNR_perTFFB(n,:)] = ...
        vicente2020_betterearsnrframe(FS, fc, [BufferedLeft(:,n) BufferedRight(:,n)], TargetLeftSpectrum, TargetRightSpectrum, InternalNoise_L, InternalNoise_R, ceiling,weightings);
end

% Compute Binaural unmasking advantage
L = length(HannWindowBU);
OverlapSample = WindowOverlap*length(HannWindowBU);
BufferedLeft = buffer(MaskerSig(:,1), L, OverlapSample, 'nodelay') .* HannWindowBU;
BufferedRight = buffer(MaskerSig(:,2), L, OverlapSample, 'nodelay') .* HannWindowBU;

weighted_BUAdv_frames = zeros(size(BufferedLeft,2),1);
BUAdv_perTFFB = zeros(size(BufferedLeft,2),length(fc));

for n = 1:size(BufferedLeft,2)
    [weighted_BUAdv_frames(n),BUAdv_perTFFB(n,:)] = ...
        vicente2020_buadvantage(FS, fc, [BufferedLeft(:,n) BufferedRight(:,n)], TargetLeftSpectrum, TargetRightSpectrum, TargetInterauralPhase, InternalNoise_L, InternalNoise_R, weightings);
end

BetterEarSNR = mean(weighted_BE_SNR_frames);
BinauralUnmaskingAdvantage = mean(weighted_BUAdv_frames);
BinauralRatio = BetterEarSNR + BinauralUnmaskingAdvantage;
end

function [left_spectrum, right_spectrum, interaural_phase ] = local_gettargetstats2(sig,fc,fs)
%Compute the (left and right) spectrum and interaural phase of the input signal sig (stereo files=2-colum matrix) sampled at fs for
%each (gammatone) frequency band with center frequencies given by fc
%Computations are similar to those used in lavandier2022.m

interaural_phase = zeros(size(fc));
left_spectrum = zeros(size(fc));
right_spectrum = zeros(size(fc));

for n = 1:length(fc)
    % filter target
    sig_left = auditoryfilterbank(sig(:,1),fs,fc(n), 'lavandier2022');    
    sig_right = auditoryfilterbank(sig(:,2),fs,fc(n), 'lavandier2022');   
    [interaural_phase(n), ~] = local_do_xcorr(sig_left,sig_right,fs,fc(n)); % cross-correlate
    % spectrum in dB based on rms of the signals (independent of signal length but not of 0 padding) rms=10*Log10(mean(sig.*sig))
    left_spectrum(n) = 10*log10(mean(sig_left.^2));
    right_spectrum(n) = 10*log10(mean(sig_right.^2));        
end
end

function [phase, coherence] = local_do_xcorr(left, right, fs, fc)
    [iacc, lags] = xcorr(left,right,round(fs/(fc*2)),'coeff'); %round(fs/(fc*2)) is for conformity with Durlach's 1972 formulation which allows time delays up to 
                                                               %+/- half the period of the channel centre frequency.
    [coherence, delay_samp] = max(iacc);
    phase = fc*2*pi*lags(delay_samp)/fs;
end

