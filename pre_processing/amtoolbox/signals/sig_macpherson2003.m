function sig = sig_macpherson2003(varargin)
%sig_macpherson2003 - Stimulus from Macpherson and Middlebrooks (2003)
%   Usage: sig = sig_macpherson2003(density,depth,phase,dur,fs,flow,fhigh,fmin)
%
%   Noise burst with sinusoidal spectral magnitude ripple.
%
%   Input parameters:
%     density : ripples/oct. Default is 1 r/o.
%     depth   : ripple depth (peak-to-trough) in dB. Default is 40 dB.
%     phase   : ripple phase in rad. Default is 0.
%     dur     : signal duration in seconds. Default is 0.25 seconds.
%     fs      : sampling rate in Hz (only required for low-pass filtering).
%               Default is 48 kHz.
%     flow    : lower corner frequency of ripple modification in Hz. Default is 1 kHz.
%     fhigh   : upper corner frequency of ripple modification in Hz. Default is 16 kHz.
%     fmin    : frequency limit of signal energy in Hz. Default is 600 Hz.
%
%   Output parameters:
%     sig     : spectrally rippled noise.
%
%   References:
%     E. A. Macpherson and J. C. Middlebrooks. Vertical-plane sound
%     localization probed with ripple-spectrum noise. jasa, 114(430), 2003.
%     
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/signals/sig_macpherson2003.php

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

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

definput.keyvals.density = 1; % ripples/oct
definput.keyvals.depth =   40; % ripple depth (peak-to-trough) in dB
definput.keyvals.phase =   0; % ripple phase in rad.
definput.keyvals.fs = 48e3; % sampling rate
definput.keyvals.fmin = 600;   % frequency limit of signal energy in Hz
definput.keyvals.flow = 1e3;   % lower corner frequency of ripple modification in Hz
definput.keyvals.fhigh = 16e3; % upper corner frequency of ripple modification in Hz
definput.keyvals.dur = 0.250; % duration

[flags,kv]  = ltfatarghelper({'density','depth','phase','dur','fs','flow','fhigh'},definput,varargin);

%% Stimulus: 
% 250-ms bursts, 20-ms raised-cosine fade in/out, flat from 0.6-16kHz
 
noiseSig = noise(kv.dur*kv.fs,1);
ph = angle(fftreal(noiseSig));

Nf = kv.dur*kv.fs/2+1;    % # positive frequency bins

f = 0:kv.fs/2/Nf:kv.fs/2;	% frequency bins
id600 = find(f<=kv.fmin,1,'last'); % index of 600 Hz (lower corner frequency of stimulus energy)
idlow = find(f<=kv.flow,1,'last'); % index of flow (ripples)
idhigh = find(f>=kv.fhigh,1,'first');  % index of fhigh (ripples)
N600low = idlow - id600 +1;   % # bins without ripple modification
Nlowhigh = idhigh - idlow +1; % # bins with ripple modification     % 
O = log2(f(idlow:idhigh)/1e3);   % freq. trafo. to achieve equal ripple density in log. freq. scale

% Raised-cosine "(i.e., cos^2)" ramp 1/8 octave wide
fup = f(idlow)*2^(1/8);       % upper corner frequency of ramp upwards 
idup = find(f<=fup,1,'last');
Nup = idup-idlow+1;
rampup = cos(-pi/2:pi/2/(Nup-1):0).^2;
fdown = f(idhigh)*2^(-1/8);  % lower corner frequency of ramp downwards
iddown = find(f>=fdown,1,'first');
Ndown = idhigh-iddown+1;
rampdown = cos(0:pi/2/(Ndown-1):pi/2).^2;
ramp = [rampup ones(1,Nlowhigh-Nup-Ndown) rampdown];
ramp = [-inf*ones(1,id600-1) zeros(1,N600low) ramp -inf*ones(1,Nf - idhigh -1)];

TF = zeros(Nf,1);  
TF(idlow:idhigh) = (kv.depth/2* sin(2*pi*kv.density*O+ kv.phase));
TF = ramp' .* TF;
TF(isnan(TF)) = -100;
sig = ifftreal(10.^(TF/20).*exp(1i*ph),2*(Nf-1));
% sexp1 = circshift(sexp1,Nf);  % IR corresponding to ripple modification

end
