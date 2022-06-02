function [dataOut] = exp_moore1997(varargin)
%EXP_MOORE1997 Figures from Moore et al. (1997)
%
%   Usage:
%     [dataOut] = exp_moore1997(flags)
%
%   The following flags can be specified:
%     'fig2'    Reproduce Fig. 2 of moore1997.
%
%     'fig3'    Reproduce Fig. 3 of moore1997.
%
%     'fig5'    Reproduce Fig. 5 of moore1997.
%
%     'fig6'    Reproduce Fig. 6 of moore1997.
%
%     'fig7'    Reproduce Fig. 7 of moore1997.
%
%   Examples:
%   ---------
%
%   To display Fig. 2 from Moore et al. (1997) use :
%
%     exp_moore1997('fig2');
%
%   To display Fig. 5 from Moore et al. (1997) use :
%
%     exp_moore1997('fig5');
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/experiments/exp_moore1997.php

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

%   AUTHOR: Thomas Deppisch
%   15.05.2020: PM: Tested, not calculated figures removed.

%% Retrieve and compute model paramters
    % Set flags

    definput.flags.type = {'missingflag','fig2','fig3','fig5','fig6','fig7'};

    [flags,kv]  = ltfatarghelper({},definput,varargin);

    if flags.do_missingflag
           flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
           sprintf('%s or %s',definput.flags.type{end-1},...
           definput.flags.type{end})];
           error('%s: You must specify one of the following flags: %s.',...
                 upper(mfilename),flagnames);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig2
        fs = 32000;
        fVec = 20:fs/2;
        data = data_glasberg2002('tfOuterMiddle1997','fieldType','free','fVec',fVec);
        figure
        semilogx(data.fOuter,data.tfOuter);
        xlim([20,20000])
        ylim([-5,20])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig3
        fs = 32000;
        fVec = 20:fs/2;
        data = data_glasberg2002('tfOuterMiddle1997','fieldType','free','fVec',fVec);
        figure
        semilogx(data.fMiddle,abs(data.tfMiddle));
        xlim([20,20000])
        ylim([0,40])
        xlabel('Frequency, kHz')
        ylabel('Middle ear effective attenuation')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig5
        % todo: figure not yet the same as in moore1997 paper
        fs = 32000;
        t = linspace(0,1,fs);
        sig = sin(2*pi*1000*t).';
        figure
        for ii=1:9
        inSig = scaletodbspl(sig,(ii+1)*10, 100);  % plot from 20 to 100 dB in 10 dB steps
        results = moore1997(inSig,fs);
        plot(results.erbN,results.eLdB,'k')
        hold on
        end
        xlim([5,40])
        ylim([0,100])
        xlabel('Number of ERB')
        ylabel('Excitation level, dB')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig6
        erbStep = 0.25;    % according to moore1997, glasberg2002
        erbFcMin = 50;
        erbFcMax = 15000;
        erbNMin = fc2erb(erbFcMin);
        erbNMax = fc2erb(erbFcMax);
        erbN = erbNMin:erbStep:erbNMax;    % numbers of erb bands
        erbFc = erb2fc(erbN); 
        dataSL = data_glasberg2002('specLoud','fVec',erbFc);
        gdB = dataSL.g;    % low level gain in cochlea amplifier
        alpha = dataSL.alpha;    % compressive exponent
        figure
        plot(gdB,alpha)
        xlim([-25,0])
        ylim([0.2,0.3])
        xlabel('G, dB')
        ylabel('alpha')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig7
        erbStep = 0.25;    % according to moore1997, glasberg2002
        erbFcMin = 50;
        erbFcMax = 15000;
        erbNMin = fc2erb(erbFcMin);
        erbNMax = fc2erb(erbFcMax);
        erbN = erbNMin:erbStep:erbNMax;    % numbers of erb bands
        erbFc = erb2fc(erbN); 
        dataSL = data_glasberg2002('specLoud','fVec',erbFc);
        gdB = dataSL.g;    % low level gain in cochlea amplifier
        a = dataSL.a;    % parameter for linearization around absolute threshold
        figure
        plot(gdB,a)
        xlim([-25,0])
        ylim([4,10])
        xlabel('G, dB')
        ylabel('A')
    end
    

end

