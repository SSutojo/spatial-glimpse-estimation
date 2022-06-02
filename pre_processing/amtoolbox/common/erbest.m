function [BW, QERB]= erbest(impulse_ch, fsig, Fs)
%erbest takes the impulse response takes its FFT and calculates its
%   ERB (Equivalent Rectangular Bandwidth) and quality factor (QERB).
%   Called by demo_lyon2011. 
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/common/erbest.php

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
le=length(impulse_ch);
F=(2*abs(fft(impulse_ch))/le).^2;
[M,Mn]=max(F(1:floor(le/2)));
E=sum(F(1:floor(le/2)));
BW=(E/M)*Fs/le;
QERB=fsig/BW;

