function [wSNR] = prudhomme2020(target,targSpec,masker,meanf0,fs, jitter,fc_target)
%PRUDHOMME2020 Compute the effective SNR taking into account harmonic cancellation
%   Usage: [predicted_SNR, BE, BU] = prudhomme2020(target,targSpec,masker,meanf0,fs, jitter,fc_target)
%
%   Input parameters:
%     target          : target
%     targSpec        : target spectrum
%     masker          : masker
%     meanf0          : mean fundamental frequency [Hz]
%     fs              : sampling frequency [Hz]
%     jitter          : jitter
%     fc_target       : center frequency of the target [Hz]
%
%   Output parameters:
%     wSNR            : weighted Signal-to-Noise Ratio
%
%   PRUDHOMME2020 monaural model computing the effective SNR taking into 
%   account harmonic cancellation, takes the target and interferer signals 
%   (sampled at fS) as input, along with masker F0 and jitter info
%
%   See also: lavandier2022 vicente2020nh vicente2020 prudhomme2020 leclere2015
%   jelfs2011
%
%   References:
%     M. Lavandier. A series of speech intelligibility models in the auditory
%     modeling toolbox. actaunited, 2022.
%     
%     P. et al. A harmonic-cancellation-based model to predict speech
%     intelligibility against a harmonic masker. jasa, 148(5):3246--3254,
%     2020.
%     
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/models/prudhomme2020.php

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
%   #Verification: Verified
%   #Requirements: MATLAB

%   AUTHOR: Matthieu Lavandier
%   adapted for AMT by Clara Hollomey (2021)


ceiling = 40;
freq_limit = 5000; %set the frequency limit to 5000 Hz
nerbs = 1:0.5:round(f2erbrate(fs/2));
fc = zeros(size(nerbs));
if length(fc)~= length(fc_target)                %check that fc for target and masker are the same
    disp('Target and masker stats should be computed at the same frequency') 
end
SNR = zeros(size(nerbs));
coeffBW = 0.6;

jittered_f0 = meanf0 + jitter*meanf0; 
    
%design the comb filter
d = fdesign.comb('notch','L,BW,GBW,Nsh',round(fs/jittered_f0),coeffBW*jittered_f0,-4,4,fs);
Hd = design(d,'SystemObject',true);
% apply the comb filter to the masker and the target
masker_HC = filter(Hd.Numerator,Hd.Denominator,masker); 
target_HC = filter(Hd.Numerator,Hd.Denominator,target);        
[targSpec_HC,~] = local_get_spectrum(target_HC,fs);

for n = 1:length(nerbs)   
    fc(n) = round(erbrate2f(nerbs(n)));         
    masker_gammatone = auditoryfilterbank(masker,fs,fc(n), 'lavandier2022'); 
    masker_HC_gammatone = auditoryfilterbank(masker_HC,fs,fc(n), 'lavandier2022'); 
    SNRwoHC = targSpec(n)-10*log10(mean(masker_gammatone.^2)); %compute the SNR without harmonic cancellation
    SNRHC = targSpec_HC(n)-10*log10(mean(masker_HC_gammatone.^2)); %compute the SNR with harmonic cancellation
    if fc(n)<freq_limit
        SNR(n) = min(ceiling,max(SNRwoHC,SNRHC)); % if the frequency is below the limit, the SNR is the best between SNR with and without harmonic cancellation
    else
        SNR(n) = min(ceiling,SNRwoHC); % if the frequency is above the limit, the SNR is the SNR without harmonic cancellation
    end
end

%integration accross frequency using SII weightings
weightings = f2siiweightings(fc);
wSNR = sum(SNR.*weightings');

end

function [spectrum, fc] = local_get_spectrum(sig,fs)

%to be used with prudhomme2020.m
%Compute the (left and right) spectrum of the input signal sig (stereo files=2-colum matrix) sampled at fs for
%each (gammatone) frequency band with center frequency given by fc
%Computations are similar to those used in lavandier2022.m

nerbs = 1:0.5:round(f2erbrate(fs/2));
fc = zeros(size(nerbs));
spectrum = zeros(length(nerbs),size(sig,2));

for n = 1:length(nerbs)
    % get filter center frequency
    fc(n) = round(erbrate2f(nerbs(n)));         
    
    for i = 1:size(sig,2)    
        % filter target
        sig_gam = auditoryfilterbank(sig(:,i),fs,fc(n), 'lavandier2022');    
        % spectrum in dB based on rms of the signals (independent of signal length but not of 0 padding) rms=10*Log10(mean(sig.*sig))
        spectrum(n,i) = 10*log10(mean(sig_gam.^2));
    end    

end
end

