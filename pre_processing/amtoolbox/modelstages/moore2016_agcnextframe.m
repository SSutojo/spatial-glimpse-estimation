function out = moore2016_agcnextframe( dLastFrame, dThisInput, aA, aR )
%MOORE2016_AGCNEXTFRAME adjusts successive short term loudness frames
%
%   Input parameters:
%     dLastFrame : last idx
%     dThisInput : current idx
%     aA : attack coefficient
%     aR : release coefficient
%
%   Output parameters:
%     out : adjusted loudness frame
%
%   version for TVL 2016 based on ANSI S3.4-2007 and Moore & Glasberg (2007)
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/modelstages/moore2016_agcnextframe.php

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

if ( dThisInput > dLastFrame )                            % attack
    out = aA * dThisInput + ( 1 - aA ) * dLastFrame;
else                                                      %release
    out = aR * dThisInput + ( 1 - aR ) * dLastFrame;
end

