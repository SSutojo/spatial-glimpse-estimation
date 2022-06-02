function exp_mckenzie2021(varargin)
%EXP_MCKENZIE2021 experiments from McKenzie et al
%
%   Usage:
%     exp_mckenzie2021('fig10d'); % (to plot Figure 10 d).
%     exp_mckenzie2021('fig10a'); % (to plot Figure 10 a).
%
%   Reproduces the plots in the listening test section (Figure 10a-d) in the 
%   Acta Acustica paper: McKenzie, T., Armstrong, C., Ward, L., Murphy, D. T., 
%   & Kearney, G. (2021). A Perceptually Motivated Spectral Difference Model 
%   for Binaural Signals. Acta Acustica (in review).
%
%   Read in perceptual listening test results and compare correlation of
%   median results to spectral difference values between the reference and
%   test stimuli. The listening test results are from the following paper:
%   McKenzie, T., Murphy, D. T., & Kearney, G. C. (2018). Diffuse-Field
%   Equalisation of Binaural Ambisonic Rendering. Applied Sciences, 8(10).
%   https://doi.org/10.3390/app8101956
%   The test compares binaural Ambisonic sounds with and without
%   diffuse-field equalisation to HRTF convolution sounds.
%
%   PEAQ and CLL are commented out in this script so that it runs without any
%   additional files necessary. To produce the respective data, download the
%   PEAQ and CLL code:
%   https://github.com/NikolajAndersson/PEAQ and
%   www.acoustics.hut.fi/-301ville/software/auditorymodel/.
%   The test sounds (ts.... etc) must be saved as wav files to work with the
%   PEAQ model.
%
%   Examples:
%   ---------
%
%   To display Figure 10a use :
%
%     exp_mckenzie2021('fig10a');
%
%   To display Figure 10d use :
%
%     exp_mckenzie2021('fig10d');
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/experiments/exp_mckenzie2021.php

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

%   Authors:
%   Thomas McKenzie, Cal Armstrong, Lauren Ward, Damian Murphy, Gavin Kearney
%   Correspondence to thomas.mckenzie@aalto.fi (happy to answer any questions
%   if you're having trouble!)

definput.flags.type={'missingflag','fig10a', 'fig10b', 'fig10c', 'fig10d'};
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

%% Read in listening test stimuli
data = amt_load('mckenzie2021', 'sig_mckenzie2021.mat');
fs = data.fs;
rsH = data.rsH;
testDirections = data.testDirections;
tsA1 = data.tsA1;
tsA3 = data.tsA3;
tsA5 = data.tsA5;
tsD1 = data.tsD1;
tsD3 = data.tsD3;
tsD5 = data.tsD5;

% combine stimuli into one matrix
ts = cat(3,tsA1,tsA3,tsA5,tsD1,tsD3,tsD5);
rs = cat(3,rsH,rsH,rsH,rsH,rsH,rsH);
tsP = permute(ts,[1 3 2]);
rsP = permute(rs,[1 3 2]);

%% Read in listening test results
data = amt_load('mckenzie2021', 'data_mckenzie2021.mat', 'data_subs');

% arrange results from repeated stimuli
resultsRaw(:,:,1) = [data.data_subs.Scen1 ; data.data_subs.Scen9];
resultsRaw(:,:,2) = [data.data_subs.Scen2 ; data.data_subs.Scen10];
resultsRaw(:,:,3) = [data.data_subs.Scen3 ; data.data_subs.Scen11];
resultsRaw(:,:,4) = [data.data_subs.Scen4 ; data.data_subs.Scen12];
resultsRaw(:,:,5) = [data.data_subs.Scen5 ; data.data_subs.Scen13];
resultsRaw(:,:,6) = [data.data_subs.Scen6 ; data.data_subs.Scen14];
resultsRaw(:,:,7) = [data.data_subs.Scen7 ; data.data_subs.Scen15];
resultsRaw(:,:,8) = [data.data_subs.Scen8 ; data.data_subs.Scen16];

% concatenate test results and order them
listTestResultsMedian = round(squeeze(median(resultsRaw))');
listTestResults = [listTestResultsMedian(:,4,:);listTestResultsMedian(:,6,:);listTestResultsMedian(:,8,:);listTestResultsMedian(:,5,:);listTestResultsMedian(:,7,:);listTestResultsMedian(:,9,:)];

%% Run spectral difference calculations
% parameters
freqRange = [20 20000]; nfft = length(rs(:,1,1));

%switch figureNumber
if flags.do_fig10a
        % BASIC SPECTRAL DIFFERENCE
        tsF = fftMatrix(tsP, fs, nfft, freqRange);
        rsF = fftMatrix(rsP, fs, nfft, freqRange);
        
        % get single values of spectral difference for all stimuli
        BSpecDiff = squeeze(mean(abs(tsF-rsF)));
        BavgSpecDiffS = mean(BSpecDiff,2);
        
elseif flags.do_fig10b
          amt_disp('perceptual evaluation of audio quality not implemented.');
%         % PERCEPTUAL EVALUATION OF AUDIO QUALITY
%         for i = 1:length(testDirectory)
%             [odg_tsA1(i), movb_tsA1(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'Ambi_1__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsA3(i), movb_tsA3(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'Ambi_3__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsA5(i), movb_tsA5(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'Ambi_5__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsD1(i), movb_tsD1(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'DFC_1__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsD3(i), movb_tsD3(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'DFC_3__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsD5(i), movb_tsD5(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'DFC_5__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%         end
%         
%         % get single values of spectral difference for all stimuli
%         QavgSpecDiffS = [odg_tsA1 odg_tsA3 odg_tsA5 odg_tsD1 odg_tsD3 odg_tsD5]';
%         
elseif flags.do_fig10c
         amt_disp('composite loudness level not implemented.');
%         % COMPOSITE LOUDNESS LEVEL
%         for i = 1:length(tsP(1,:,1))
%             CLL_difference(i,:) = CLL(tsP(:,i,:),rsP(:,i,:),fs);
%         end
%         % get single values of spectral difference for all stimuli
%         CavgSpecDiffS = mean(abs(CLL_difference),2);
%         
elseif flags.do_fig10d
        % PERCEPTUAL SPECTRAL DIFFERENCE
        f.fs = fs;f.nfft = nfft;f.minFreq = freqRange(1); f.maxFreq = freqRange(2);
        [~,PSpecDiff] = mckenzie2021(tsP,rsP,0,f,0); %no fft pre model
        PSpecDiff = squeeze(PSpecDiff);
        
        % get single values of spectral difference for all stimuli
        PavgSpecDiffS = mean(PSpecDiff,2);
end

%% Calculate Pearson's Correlation Coefficient
% between spectral difference values and perceptual listening test results
%switch figureNumber
if flags.do_fig10a
        [rbsd, pbsd] = corrcoef(BavgSpecDiffS,listTestResults);
        disp(strcat('BSD correlation=',num2str(rbsd(2,1)),', p=',num2str(pbsd(2,1))));
        
%     case 'fig10b'
%         [rqsd, pqsd] = corrcoef(QavgSpecDiffS,listTestResults);
%         disp(strcat('PEAQ correlation = ',num2str(rqsd(2,1)),', p = ',num2str(pqsd(2,1))));
%         
%     case 'fig10c'
%         [rcsd, pcsd] = corrcoef(CavgSpecDiffS,listTestResults);
%         disp(strcat('CLL correlation = ',num2str(rcsd(2,1)),', p = ',num2str(pcsd(2,1))));
        
elseif flags.do_fig10d
        [rpsd, ppsd] = corrcoef(PavgSpecDiffS,listTestResults);
        disp(strcat('PSD correlation=',num2str(rpsd(2,1)),', p=',num2str(ppsd(2,1))));
end

%% Plot PCC vs Listening Test Results
%switch figureNumber
if flags.do_fig10a
        h = figure; scatter(listTestResults,BavgSpecDiffS,60,'x','LineWidth',1.5,'MarkerEdgeColor','k');
        ylabel('BSD (dB)'); xlabel('Median Results');
        set(gca,'FontSize', 14); set(gcf, 'Color', 'w');
        pbaspect([1.7 1 1]); grid on; box on;
        % saveFig(h,'bsd_dfe.pdf',2)
        
%     case 'fig10b'
%         h = figure; scatter(listTestResults,QavgSpecDiffS,60,'x','LineWidth',1.5,'MarkerEdgeColor','k');
%         ylabel('PEAQ (ODG)'); xlabel('Median Results');
%         set(gca,'FontSize', 14); set(gcf, 'Color', 'w');
%         pbaspect([1.7 1 1]); grid on; box on;
%         % saveFig(h,'qsd_dfe.pdf',2)
%         
%     case 'fig10c'
%         h = figure; scatter(listTestResults,CavgSpecDiffS,60,'x','LineWidth',1.5,'MarkerEdgeColor','k');
%         ylabel('CLL (dB)'); xlabel('Median Results');
%         set(gca,'FontSize', 14); set(gcf, 'Color', 'w');
%         pbaspect([1.7 1 1]); grid on; box on;
%         % saveFig(h,'csd_dfe.pdf',2)
        
elseif flags.do_fig10d
        h = figure; scatter(listTestResults,PavgSpecDiffS,70,'x','LineWidth',1.5,'MarkerEdgeColor','k');
        ylabel('PSD (sones)'); xlabel('Median Results');
        set(gca,'FontSize', 14); set(gcf, 'Color', 'w');
        pbaspect([1.7 1 1]); grid on; box on;
        % saveFig(h,'psd_dfe.pdf',2)
end

%% Extra functions
    function [matrix_output_fft, freq_vector_fft,fft_abs_matrix_input] = fftMatrix(matrix_input, Fs, Nfft, freq_range)
        % Function to calculate the single sided frequency spectrum of two matrices
        % of HRIRs for a specified frequency range. Returns FFT of input matrix as
        % the absolute FFT in dB for the specified frequency range with the
        % associated frequency vector.
        
        % Take FFT of matrices
        fft_matrix_input = fft(matrix_input, Nfft); % Get Fast Fourier transform
        
        % Compute freq bins for x-axis limits
        fr_low = round(freq_range(1)*Nfft/Fs);
        fr_high = round(freq_range(2)*Nfft/Fs);
        
        % Get absolute values for frequency bins
        fft_abs_matrix_input = abs(fft_matrix_input(fr_low:fr_high,:,:));
        
        % Get values in dB
        matrix_output_fft = 20*log10(fft_abs_matrix_input);
        
        % Frequency vector for plotting
        f = 0:Fs/Nfft:Fs-(Fs/Nfft);
        freq_vector_fft = f(fr_low:fr_high);
    end

% Composite loudness level difference
    function [CLL_difference,freqs] = CLL(input,ref,Fs)
        if length(input(:,1)) < length(input(1,:))
            input = input';
        end
        if length(ref(:,1)) < length(ref(1,:))
            ref = ref';
        end
        
        % requires Karjalainen Auditory Toolbox.
        [~, ~, ~, CLL_input, freqs] = simuspatcues_KarAudMod(input(:,1),input(:,2),Fs);
        [~, ~, ~, CLL_ref] = simuspatcues_KarAudMod(ref(:,1),ref(:,2),Fs);
        
        CLL_difference = CLL_input-CLL_ref;
    end

end

