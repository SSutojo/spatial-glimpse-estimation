function amt_stop
%amt_stop Removes all amt_related paths and clears persistent variables
%
%   'amt_stop' removes all paths containing the 'amt_basepath' and clears the
%   memory of all persistent variables within the amt core functions. This 
%   function facilitates working with different AMT versions.
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/core/amt_stop.php

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

%   #Author: Clara Hollomey (2021)

%retrieve all current paths
allpaths = path;
S = regexp(allpaths, pathsep, 'split');

if exist('arg_amt_configuration', 'file')
    %find the AMT topfolder
    [flags, kv] = amt_configuration;
    if isempty(flags)
        disp('No configuration available. AMT may already have been uninstalled.')
        return;
    end
    %pathelements = regexp(kv.path, filesep, 'split');
    %topfolder = pathelements{numel(pathelements)};
    topfolder = kv.path;
    topfolder = topfolder(1:end-1);
    if isempty(topfolder)
        %just make sure to really capture the name of the topfolder (system-dependent)
        topfolder = pathelements{numel(pathelements)-1};    
    end

    %clear persistent variables from all core functions
    functions = dir([kv.path, filesep, 'core', filesep, 'amt_*']);

    for ii = 1:numel(functions)
      try
      munlock(functions(ii).name);
      catch
      end
      clear(char(functions(ii).name));
    end
%    amt_configuration('amtrunning', 0);
    %delete all the paths that contain the topfolder name
    deletepaths = strfind(S, topfolder);
    deletelines=find(cellfun(@(c) ~isempty(c), deletepaths));

    if sum(deletelines) > 0
        rmpath(S{deletelines});
        %amt_disp not possible at this stage anymore as path may already
        %have been removed
        if flags.do_verbose, disp(['All AMT ' kv.version{1}(5:end) ' paths were removed.']); end
    else
        if flags.do_verbose, disp('No paths to remove.'); end
    end
else
    disp('AMT not loaded.')
end

