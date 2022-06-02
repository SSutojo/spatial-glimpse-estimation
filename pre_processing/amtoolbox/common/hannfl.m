function h = hannfl(len,h1len,h2len);
%HANNFL plots a hann window 
%   the symmetry of the window can be controlled by h1len and h2len
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/common/hannfl.php

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

if len > 0
h = ones(len,1);

switch h1len
case 0
otherwise
   h(1:h1len)=(1-cos(pi/h1len*[0:h1len-1]))/2;
end

switch h2len
case 0
otherwise
   h(end-h2len+1:end)=(1+cos(pi/h2len*[1:h2len]))/2;
end

else
end

