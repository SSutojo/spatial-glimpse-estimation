function nlf = lyon2011_ohcnlf(velocities, CAR_coeffs)
%lyon2011_ohcnlf velocity after OHC processing
%
%   Usage:
%     nlf = lyon2011_ohcnlf(velocities, CAR_coeffs)
%
%   Input parameter:
%     velocities : input velocities
%     CAR        : struct, CAR coefficients
%
%   Output parameter:
%     nlf        : velocity after OHC processing
%
%   start with a quadratic nonlinear function, and limit it via a
%   rational function; make the result go to zero at high
%   absolute velocities, so it will do nothing there.
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/modelstages/lyon2011_ohcnlf.php

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

%   #Author: Amin Saremi (2016) adaptations for the AMT (based on <https://github.com/google/carfac>, Richard F. Lyon)
%   #Author: Clara Hollomey (2021) adaptation for the AMT 1.0
%   #License: gpl3


nlf = 1 ./ (1 + ...
  (velocities * CAR_coeffs.velocity_scale + CAR_coeffs.v_offset) .^ 2 );


