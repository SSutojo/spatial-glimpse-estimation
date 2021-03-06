function erb = fc2erb(fc)
%fc2erb calculates the ERB index from a given frequency
%
%   Usage: erb = fc2erb(fc)
%
%   converts from [Hz] to [cam]
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/common/fc2erb.php

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

    erb = 21.366*log10(0.004368.*fc + 1);
end

