%% test installation mode
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/test_amt_start.php

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
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: All messages displayed, install toolboxes');
amt_start('install');
disp('TEST AMT_START: DEMO_BAUMGARTNER2013 should run and display results now'); 
demo_baumgartner2013;
disp('TEST AMT_START: AMT_STOP should display a line'); 
amt_stop;

%% test default verbose mode
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: All messages displayed');
amt_start;
disp('TEST AMT_START: DEMO_BAUMGARTNER2013 should display results now'); 
demo_baumgartner2013;
disp('TEST AMT_START: AMT_STOP should display a line'); 
amt_stop;

%% test silent mode
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: SILENT MODE: no messages displayed');
amt_start('silent');
disp('DTEST AMT_START: DEMO_BAUMGARTNER2013 should display nothing'); 
demo_baumgartner2013;
disp('TEST AMT_START: AMT_STOP should display nothing'); 
amt_stop;

%% test documentation mode
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: DOCUMENTATION MODE should display as normal');
amt_start('documentation');
disp('TEST AMT_START: DEMO_BAUMGARTNER2013 should display nothing'); 
demo_baumgartner2013;
disp('TEST AMT_START: AMT_STOP should display nothing'); 
amt_stop;

%% test done
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: FINISHED');

