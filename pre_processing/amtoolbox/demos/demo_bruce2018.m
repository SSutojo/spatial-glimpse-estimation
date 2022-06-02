%DEMO_BRUCE2018 generate an example neurogram
%
%   Figure 1: Mean-rate neurogram
%
%   Figure 2: Fine-timing neurogram
%
%   Figure 3: S_out neurogram
%
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/demos/demo_bruce2018.php

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

% Set audiogram to simulate normal hearing
ag_fs = [125 250 500 1e3 2e3 4e3 8e3];
ag_dbloss = [0 0 0 0 0 0 0]; 
% Set audiogram for example high-freq hearing loss
% ag_fs = [125 250 500 1e3 2e3 4e3 8e3];
% ag_dbloss = [0 0 0 20 40 60 80]; 

species = 'human'; % Human cochlear tuning (Shera et al., 2002)
numL = 10; numM = 10; numH = 30;
[stim, Fs_stim] = amt_load('bruce2018','defineit.wav');
stimdb = 60; % SPL of speech (in dB)
stim = scaletodbspl(stim,stimdb);

numCF = 40;
flow = 250;
fhigh = 16e3;

fc = logspace(log10(flow), log10(fhigh),numCF);

out = bruce2018(stim, Fs_stim, fc, 'ag_fs', ag_fs, 'ag_dbloss', ag_dbloss, ...
  'numL', numL,'numM', numM,'numH', numH,'nrep',1,'outputPerCF');
amt_disp();
neurogram_ft = out.neurogram_ft;
neurogram_mr = out.neurogram_mr;
neurogram_Sout = out.neurogram_Sout;
t_ft = out.t_ft;
t_mr = out.t_mr;
t_Sout = out.t_Sout;
CFs = out.fc;

ng1=figure;
set(ng1,'renderer','painters');
winlen = 256; % Window length for the spectrogram analyses
sp1 = subplot(2,1,1);
[s,f,t] = specgram([stim; eps*ones(round(t_mr(end)*Fs_stim)-length(stim),1)],winlen,Fs_stim,winlen,0.25*winlen);
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
xlabel('Time');
ylabel('Frequency (kHz)');
title('Spectrogram')
xl = xlim;
sp2=subplot(2,1,2);

plot_bruce2018(t_mr,CFs,neurogram_mr,sp2);

caxis([0 80])
title('Mean-rate Neurogram')
xlim(xl)

ng2=figure;
set(ng2,'renderer','painters');
winlen = 256; % Window length for the spectrogram analyses
sp1 = subplot(2,1,1);
[s,f,t] = specgram([stim; eps*ones(round(t_mr(end)*Fs_stim)-length(stim),1)],winlen,Fs_stim,winlen,0.25*winlen);
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
xlabel('Time');
ylabel('Frequency (kHz)');
title('Spectrogram')
xl = xlim;
sp2=subplot(2,1,2);
plot_bruce2018(t_ft,CFs,neurogram_ft,sp2);
caxis([0 20])
title('Fine-timing Neurogram')
xlim(xl)

ng3=figure;
set(ng3,'renderer','painters');
winlen = 256; % Window length for the spectrogram analyses
sp1 = subplot(2,1,1);
[s,f,t] = specgram([stim; eps*ones(round(t_mr(end)*Fs_stim)-length(stim),1)],winlen,Fs_stim,winlen,0.25*winlen);
imagesc(t,f/1e3,20*log10(abs(s)/sum(hanning(winlen))*sqrt(2)/20e-6));
axis xy; axis tight;
hcb = colorbar;
set(get(hcb,'ylabel'),'string','SPL')
caxis([stimdb-80 stimdb])
ylim([0 min([max(CFs/1e3) Fs_stim/2e3])])
xlabel('Time');
ylabel('Frequency (kHz)');
title('Spectrogram')
xl = xlim;
sp2=subplot(2,1,2);
plot_bruce2018(t_Sout,CFs,neurogram_Sout*diff(t_Sout(1:2)),sp2);
caxis([0 6])
title('S_{out} Neurogram')
xlim(xl)

