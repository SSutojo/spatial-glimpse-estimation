%DEMO_RELANOIBORRA2019 experiments from Relano-Iborra et al. 2019
%
%   Usage: demo_relanoiborra2019;
%
%   This script reproduces Figure 3 of Relano-Iborra et al. It finds the free 
%   parameters of the logistic function that relates the output of the sCASP 
%   model (correlation coefficient d) with the percentage of correct answers. 
%   The data used for this fitting is taken from Nielsen and Dau (2009). The 
%   testing conditions are replicated and fed to the sCASP model. Later a least 
%   squares analysis is used to fit the function to the data.
%
%   figure
%
%      sCASP mapping for the CLUE speech corpus
%
%   AUTHOR: Helia Relano Iborra, September 2020.
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/demos/demo_relanoiborra2019.php

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


%% Initialization

Pcorrect_human = [0 8 35 71 90 100 ]; %% Human data
SNRs= -8:2:2;

x = amt_load('relanoiborra2019', 'single_150_SentArray22kHz_varLength.mat');
sentenceArray = x.sentenceArray;
fsSent = 22050;

noise_name = 'SSN_CLUE_22kHz.wav'; % SSN Noise
speechSPL = 65;
Nsentences=5;

flow = 100;
fhigh = 8e3;        
sbj = 'NH';
fsRef = 22050; % sampling freq.
Pref= 20e-6; % Transformation to Pascals

d= zeros(length(SNRs), Nsentences);

%% Run experiment:

for q=1:Nsentences
    
    amt_disp(['Processing sentence: ' num2str(q) ' out of ' num2str(Nsentences)]);
    speech  = sentenceArray{q};
    speech = resample(speech, fsRef, fsSent );
    
    N_org= length(speech);  % Calculate length of original sentence
    speech = [speech; speech]; % Prepane the same sentence
    
    speech = Pref*speech*(1/rms(speech))*10^((speechSPL)/20); % Set speech level
    N = length(speech);  % Overall speech size
    
    
    for n=1:length(SNRs)
        
        %noise = audioread(noise_name);
        noise = amt_load('relanoiborra2019', noise_name);
        Nsegments = floor(length(noise)/N);
        startIdx = randi(Nsegments-2 ,1)*N;
        noise = noise(startIdx:startIdx+N -1)'; % random segment from the noise file
        
        noise = Pref*(noise./rms(noise)*10^((speechSPL-SNRs(n))/20)); % sets the level of the noise signal
        
        if size(noise) ~= size(speech)
            noise = noise'; % Flips noise signal if needed
        end
        
        test = noise + speech; % Gerating the noisy mixture               

        tmp = relanoiborra2019(speech, test, fsRef, flow, fhigh, 'N_org', N_org, 'subject', sbj);
        
        d(n, q) = tmp.dfinal;  % correlation value per sentence and SNR
        
    end % End loop over SNRs
end % End loop over sentences

d_mean= mean(d, 2); % Model outputs averaged across sentences

%% Fitting

xdata = d_mean';
ydata = Pcorrect_human;

fun = @(a,xdata) 100./(1 + exp(a(1)*xdata + a(2))); %Logistic function

pguess = [0 0]; %starting guess
[fit_param,R,J,CovB,MSE] = nlinfit(xdata,ydata,fun,pguess); % non-linear least squares optimization

% Goodness of fit
ysim= fun(fit_param, xdata);
rsq2 = 1 - sum(R.^2) / sum((ydata- mean(ydata)).^2);
varagout{1} = 0;
%% Plotting

x = linspace(0, 1, 200);
fit_funct= fun(fit_param, x);

figure
scatter(d_mean, Pcorrect_human, 68, 'filled', 'r')
hold on
plot(x, fit_funct,'k', 'LineWidth', 2)
xlabel('Model output (d)'), ylabel('% correct'), title('sCASP mapping for CLUE material')
le= legend(' Nielsen & Dau (2009) data', 'f_C_L_U_E', 'Location', 'southeast');
set(le, 'box', 'off')
text(3,12,{[' R^2 = ',num2str(rsq2,2)]},'fontsize',12,'FontName', 'cambria');
set(gca, 'fontsize',12)

