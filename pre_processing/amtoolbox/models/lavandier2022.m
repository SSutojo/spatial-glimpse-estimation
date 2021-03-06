function [twoears_benefit, weighted_bmld, weighted_better_ear] = lavandier2022(target_in,int_in,fs)
%LAVANDIER2022 Compute the binaural 'effective' target-to-interferer ratio
%   Usage: [twoears_benefit, weighted_bmld, weighted_better_ear] = lavandier2022(target_in,int_in,fs)
%
%   Input parameters:
%     target_in       : target
%     int_in          : interferer
%     fs              : sampling frequency [Hz]
%
%   Output parameters:
%     twoears_benefit       : effective target to interferer ratio
%     weighted_bmld         : weighted binaural masking level difference
%     weighted_better_ear   : weighted better ear advantage
%
%   LAVANDIER2022 computes the binaural 'effective' target-to-interferer ratio 
%   target_in and int_in are signals produced at the ears: stereo files 
%   (2-column matrices) of the same sampling frequency fs
%
%   See also: lavandier2022 vicente2020nh vicente2020 prudhomme2020 leclere2015
%   jelfs2011
%
%   References:
%     M. Lavandier. A series of speech intelligibility models in the auditory
%     modeling toolbox. actaunited, 2022.
%     
%     L. et al. Binaural prediction of speech intelligibility in reverberant
%     rooms with multiple noise sources. jasa, 131(1):218--231, 2012.
%     
%     L. et al. Speech segregation in rooms: Monaural, binaural and
%     interacting effects of reverberation on target and interferer. jasa,
%     123(4):2237--2248, 2008.
%     
%     S. Jelfs, J. Culling, and M. Lavandier. Revision and validation of a
%     binaural model for speech intelligibility in noise. Hearing Research,
%     2011.
%     
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/models/lavandier2022.php

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


nerbs = 1:0.5:round(f2erbrate(fs/2));
fc = zeros(size(nerbs));
bmld_prediction = zeros(size(nerbs));
better_ear_prediction = zeros(size(nerbs));

for n = 1:length(nerbs)
    % get filter cf
    fc(n) = round(erbrate2f(nerbs(n)));    
    % filter target and interferer separately
    targ_left = auditoryfilterbank(target_in(:,1),fs,fc(n), 'lavandier2022');    
    targ_right = auditoryfilterbank(target_in(:,2),fs,fc(n), 'lavandier2022');   
    int_left = auditoryfilterbank(int_in(:,1),fs,fc(n), 'lavandier2022');       
    int_right = auditoryfilterbank(int_in(:,2),fs,fc(n), 'lavandier2022');  
    % BMLD
    [int_phase, int_coherence] = local_do_xcorr(int_left,int_right,fs,fc(n)); % cross-correlate
    [target_phase] = local_do_xcorr(targ_left,targ_right,fs,fc(n));    
    bmld_prediction(n) = bmld(int_coherence,target_phase,int_phase,fc(n));    
    % better-ear SNR in dB based on rms of the signals (independent of
    % signal length but not of 0 padding) rms=10*Log10(mean(sig.*sig))
    left_SNR = 10*log10(mean(targ_left.^2)/mean(int_left.^2));
    right_SNR = 10*log10(mean(targ_right.^2)/mean(int_right.^2));
    better_ear_prediction(n) = max(left_SNR,right_SNR);   
end

%integration accross frequency using SII weightings
weightings = f2siiweightings(fc);
weighted_bmld = sum(bmld_prediction.*weightings');
weighted_better_ear = sum(better_ear_prediction.*weightings');

twoears_benefit = weighted_better_ear + weighted_bmld;

end

function [phase, coherence] = local_do_xcorr(left, right, fs, fc)
    [iacc, lags] = xcorr(left,right,round(fs/(fc*2)),'coeff'); %round(fs/(fc*2)) is for conformity with Durlach's 1972 formulation which allows time delays up to 
                                                               %+/- half the period of the channel centre frequency.
    [coherence, delay_samp] = max(iacc);
    phase = fc*2*pi*lags(delay_samp)/fs;
end


