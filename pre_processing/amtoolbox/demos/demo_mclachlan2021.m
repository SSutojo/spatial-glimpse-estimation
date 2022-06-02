%DEMO_MCLACHLAN2021 Demo for dynamic full sphere localization model based on Reijniers et al. (2014)
%
%   DEMO_MCLACHLAN2021(flag) demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   on the full sphere using the localization model from Reijniers et al. (2014).
%
%   Figure 1: Baseline prediction: averaged polar and lateral accuracy
%
% 
%   This demo computes the baseline prediction (localizing broadband 
%   sounds with own ears) for an exemplary listener (NH12).
%
%      
%
%   See also: mclachlan2021 mclachlan2021 plot_reijniers2014
%               reijniers2014_metrics
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/demos/demo_mclachlan2021.php

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

% AUTHOR : Glen McLachlan

%% Settings
%num_dirs = 500; 
%[dirs,~,~,~] = ParticleSampleSphere('N',num_dirs); 
%dirs = load('dirs.mat'); dirs = dirs.dirs;
%[az, el] = cart2sph(dirs(:,1),dirs(:,2),dirs(:,3));

az = [];          % azimuth target angle in degrees
el = [];          % elevation target angle in degrees
rot_type = 'yaw'; % rotation type ('yaw','pitch','roll')
rot_size = 10;    % rotation amount in degrees
num_exp = 10;     % # of virtual experimental runs
sig_itdi = 0.8;   % variance on subsequent looks

assert(numel(az)==numel(el))

amt_disp('Experiment conditions:','documentation');
for i=1:length(az)
    amt_disp(sprintf('  Azimuth: %0.1f deg, Elevation: %0.1f deg', ...
        az(i), el(i)),'documentation');
end
amt_disp(sprintf('  Repetitions: %i', num_exp));
amt_disp('------------------------');


%% Get listener's data
SOFA_obj = SOFAload(fullfile(SOFAdbPath,'baumgartner2013', ...
                        'ARI_NH12_hrtf_M_dtf 256.sofa'));
%SOFA_obj = SOFAload('HRIR_L2354.sofa');


%% Preprocessing source information for both directions
[template, target] = mclachlan2021_preproc(SOFA_obj, ... 
                    'targ_az', az, 'targ_el', el, 'rot_type', rot_type);

%% Run virtual experiments
[doa, params] = mclachlan2021(template, target, 'num_exp', num_exp, ...
    'rot_type', rot_type, 'rot_size', rot_size, 'sig_itdi', sig_itdi);

%% Calcualte performance measures 
amt_disp('------------------------')
amt_disp('Performance Predictions:')
amt_disp('------------------------')

met = mclachlan2021_metrics(doa);
met.entropy = mean(params.entropy);
met.information = mean(params.information);

%amt_disp(sprintf('Lateral accuracy: %0.2fdeg', met.accL))
amt_disp(sprintf('Lateral RMS: %0.2f', met.rmsL))
%amt_disp(sprintf('Elevation accuracy: %0.2fdeg', met.accE))
amt_disp(sprintf('Polar RMS: %0.2f', met.rmsP))
%amt_disp(sprintf('Quadrant error: %0.6f%%',met.querr))
amt_disp(sprintf('Mean entropy: %0.2fbits', met.entropy))
amt_disp('------------------------')

met.nexp = num_exp;
met.rot_type = rot_type;
met.rot_size = rot_size;
met.sig = 1;
met.coords = target.coords;
met.post_prob = params.post_prob;

% this only works if all directions are tested
%plot_reijniers2014(target.coords,params.information);
%title('Information for yaw=0.1 in bits')

