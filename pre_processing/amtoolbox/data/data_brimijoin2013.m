function data = data_brimijoin2013
%DATA_BRIMIJOIN2013 - Data from Brimijoin et al. (PLoS One, 2013)
%   Usage: data = data_brimijoin2013
%
%   Mean externalization likelihood of NH listeners extracted from Fig. 6 
%
%   Output parameters:
%     data    : structure with fields
%                 response ... mean externalization likelihood
%                 legend ... conditions (columns of response matrix)
%                 mix ... head-present re head-absent mix
%                 xlabel ... abscissa label of original figure
%                 ylabel ... ordinate label of original figure
%
%   References:
%     W. O. Brimijoin, A. W. Boyd, and M. Akeroyd. The contribution of head
%     movement to the externalization and internalization of sounds. PLOS
%     ONE, 2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/data/data_brimijoin2013.php

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

data.mix = 0:.2:1;
data.xlabel = 'Head-present / Head-absent Mix';
data.ylabel = 'Proportion externalized responses';
data.legend = {'Head: fixed; Signal: fixed';'Head: fixed; Signal: move';...
  'Head: move; Signal: fixed';'Head: move; Signal: move'};
data.response = ...
  [.08,.09,.18,.55,.70,.78;...
   .07,.08,.29,.51,.79,.73;...
   .13,.07,.08,.17,.32,.38;...
   .27,.51,.54,.81,.89,.89]';


end
