function data = data_best2005
%DATA_BEST2005  Listener averages of absolute polar angle error and SCC
%   Usage: data = data_best2005
%
%   Output parameters:
%     data.qe     : quadrant error rates
%     data.ale    : absolute lateral angle error
%     data.ape    : absolute polar angle error
%     data.seape  : standard error of absolute polar angle error
%     data.meta   : condition labels corresponding to data entries
%
%   DATA_BEST2005 returns listener averages of absolute lateral and polar
%   angle error and SCCs from Best et al. (2005): Fig.5, Fig.10, and Tab.2
%
%
%   References best2005highfrequency
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/data/data_best2005.php

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

% AUTHOR: Robert Baumgartner, Roberto Barumerli

% Mean quadrant error rates from Tab.2 (Exp. I) and Tab. 4 (Exp. II)
data.qe =     [3.0,mean([8.4,6.7]),12.0,13.8,16.4];

% Mean absolute polar angle errors from Fig.5b (Exp. I) and Fig.10b (Exp. II)
data.ape =    [21,30,38,42.5,46];
data.seape =  ones(1,5);

% Mean absolute lateral angle errors from Fig.5a (Exp. I) and Fig.10a (Exp. II)
data.ale =    [7.87, 9.63, 9.70, 9.67, 8.60];

data.meta =   {'BB noise','BB speech','-20dB','-40dB','-60dB'};
end
