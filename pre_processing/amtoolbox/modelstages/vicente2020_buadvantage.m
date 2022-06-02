function [weighted_BUAdv,BUAdv_perFB] = vicente2020_buadvantage(fs, fc, MaskerSig, TargetLeftSpectrum, TargetRightSpectrum, TargetIPD, InternalNoise_L, InternalNoise_R, weightings)
%vicente2020_buadvantage calculates the bmld advantage
%   Usage: [weighted_BUAdv,BUAdv_perFB] = vicente2020_buadvantage(fs, fc, MaskerSig, TargetLeftSpectrum, TargetRightSpectrum, TargetIPD, InternalNoise_L, InternalNoise_R, weightings)
%
%   compute IN ALL frequency bands centered at fc, the BE SNR, after considering the audiogram, applying ceiling.
%   int_in is ear signals (in a time-frame)
%   t_left/right are (longterm) power of the target at these ears in this band, t_phase is its corresponding interaural phase
%   fs=sampling frequency
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/modelstages/vicente2020_buadvantage.php

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

BUAdv_perFB = zeros(1,length(fc));

for n = 1:length(fc)
    LeftMasker = auditoryfilterbank(MaskerSig(:,1),fs,fc(n), 'lavandier2022');        
    RightMasker = auditoryfilterbank(MaskerSig(:,2),fs,fc(n), 'lavandier2022');
    if 20*log10(rms(LeftMasker)) > InternalNoise_L(n) &&  20*log10(rms(RightMasker)) > InternalNoise_R(n) &&  TargetLeftSpectrum(n) > InternalNoise_L(n) &&  TargetRightSpectrum(n) > InternalNoise_R(n)
        [MaskerIPD, MaskerCoherence] = local_do_xcorr(LeftMasker, RightMasker, fs, fc(n)); 
        BUAdv_perFB(n) = bmld(MaskerCoherence, TargetIPD(n), MaskerIPD, fc(n));                                       
    end
end
%SII weightings+integreation across frequency
weighted_BUAdv = sum(BUAdv_perFB .* weightings');

end

function [phase, coherence] = local_do_xcorr(left, right, fs, fc)
    [iacc, lags] = xcorr(left,right,round(fs/(fc*2)),'coeff'); %round(fs/(fc*2)) is for conformity with Durlach's 1972 formulation which allows time delays up to 
                                                               %+/- half the period of the channel centre frequency.
    [coherence, delay_samp] = max(iacc);
    phase = fc*2*pi*lags(delay_samp)/fs;
end
