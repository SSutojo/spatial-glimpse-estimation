function plot_mckenzie2021(spectralDifference,testDirections)
%PLOT_MCKENZIE2021 calculate max and minimum values
%
%   Usage:
%     plot_mckenzie2021(spectralDifference,testDirections);
%
%   Input parameters:
%     spectralDifference   : vector of spectral differences
%     testDirections       : vector of test directions
%
%   This function generates a spherical heatmap.
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/plot/plot_mckenzie2021.php

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

maxValue = max(spectralDifference);
minValue = min(spectralDifference);

figure('position',[100,100,1200 400]);
for i = 1:length(spectralDifference)/length(testDirections)
    subplot(2,3,i)
    
    % plot values on spherical map
    local_heatmap_plot(testDirections(:,1),testDirections(:,2),spectralDifference((i-1)*length(testDirections)+1:i*(length(testDirections))));
    
    % title for each subplot. 1OA refers to first order Ambisonics. DFE refers to diffuse-field equalisation.
    averagePSDvalue = mean(spectralDifference((i-1)*length(testDirections)+1:i*(length(testDirections))));
    switch i
        case 1
            title(strcat('1OA no DFE. Mean PSD=',num2str(averagePSDvalue),' sones'))
        case 2
            title(strcat('3OA no DFE. Mean PSD=',num2str(averagePSDvalue),' sones'))
        case 3
            title(strcat('5OA no DFE. Mean PSD=',num2str(averagePSDvalue),' sones'))
        case 4
            title(strcat('1OA with DFE. Mean PSD=',num2str(averagePSDvalue),' sones'))
        case 5
            title(strcat('3OA with DFE. Mean PSD=',num2str(averagePSDvalue),' sones'))
        case 6
            title(strcat('5OA with DFE. Mean PSD=',num2str(averagePSDvalue),' sones'))
    end
    c2 = colorbar; c2.Label.String = 'PSD (sones)';
    caxis([minValue,maxValue]);
    xlim([0 180]); % as this data is only in a hemisphere, constrain axis limits
end

end

function local_heatmap_plot(az,el,psd)
% Plot spectral difference of a large spherical set of points on a
% rectangular plot. Need input vectors of azimuth, elevation and the
% spectral difference of the points.


xlin = linspace(min(az),max(az),180*2);
ylin = linspace(min(el),max(el),90*2);
[X,Y] = meshgrid(xlin,ylin);

Z = griddata(az,el,psd,X,Y,'cubic');

surf(X,Y,Z,'EdgeColor','none')
xlabel('Azimuth (°)'); ylabel('Elevation (°)');
set(gca, 'XDir', 'reverse', 'YTick', -75:75:75, 'XTick', -150:75:150);
xlim([-180 180]); ylim([-90 90]); view ([0 90]);
colormap(flipud(parula)); set(gcf, 'Color', 'w');
axis tight; box on; pbaspect([2 1 1]);
end
