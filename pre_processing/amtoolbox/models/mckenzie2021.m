function [diff, avgDiff_ERB, norm] = mckenzie2021( A, B, domFlag, f, norm, w, plotFlag, lim)
%MCKENZIE2021 Binaural perceptual similarity
%
%   Input parameters:
%     A       : Reference input data matrix
%     B       : Test input data matrix
%     domFlag : Specifies whether input data is in time (0) / frequency (1) / frequency_dB (2) domain
%     f       : struct, if (domFlag == 0) contains fs, nfft, minFreq, maxFreq; if (domFlag == 1 or 2) Specifies the FFT sample frequencies
%
%   Optional input parameters:
%
%     'norm'             if empty (DEFAULT), Iterate to find optimal normalisation          
%                      else, apply norm dB or normalisation to input B 
%
%     'w'                (DEFAULT=1) Sample point weightings
%
%     'plotFlag'         (DEFAULT=0) Don't show (0) or show (1) plot normalisation curve
%
%     'lim'              (DEFAULT=0.05) if lim >= 1 it specififies the 
%                      number of normalisation iterations. If lim < 1 it 
%                      specifies the resolution/change in PSD that must be    
%                      reached by sequential normalisation steps i.e. 
%                      iterations will stop once they result in changes 
%                      to the PSD value that are less than lim (or 100000
%                      iterations)
%
%     'SPL'            (DEFAULT=75) the average dB SPL value at which comparisons are made
%
%     'initInc'          (DEFAULT=0.2) The starting offset increment
%
%     'scale'            (DEFAULT=0.4) The absolute scaling value by which   
%                      the offset increment is adjusted each time an 
%                      increase in PSD in detected
%
%
%   This function compares the spectra of A and B in terms of PERCEPTUALLY
%   WEIGHTED error. For multi-dimensional inputs the comparison is made along
%   the first dimension. Averages are output for each column of data.
%   Perceptual error takes into account the lesser importance of quieter sounds
%   and less sensitive frequencies. Frequency bins are weighted with respect to
%   the ISO 226 loudness curves for an average listening level of 75 dB SPL by
%   default. A contribution half as loud is deemed half as important using the
%   Sones scale.
%   The perceptual average difference is further weighted with respect to ERB
%   bandwidth. Simply, this reduces the contribution to the average calculation
%   of higher frequency components where our ears are less sensitive. It's like
%   a logarithmic type average.
%   An iterative optimisation process is used to find the input normalisation
%   which results in the lowest error metric. This is generally somewhere
%   around the point at which the two input signals have the same mean value -
%   but can easily vary by a few dB. The optimum normalisation is different for
%   perceptual / absolute error metrics.
%
%
%   Authors:
%   Thomas McKenzie, Cal Armstrong, Lauren Ward, Damian Murphy, Gavin Kearney
%   Correspondence to thomas.mckenzie@aalto.fi (happy to answer any questions
%   if you're having trouble!)
%
%
%   See also: plot_mckenzie2021 demo_mckenzie2021 exp_mckenzie2021
%
%   References:
%     T. McKenzie, C. Armstrong, L. Ward, D. T. Murphy, and G. Kearney. A
%     perceptually motivated spectral difference model for binaural signals.
%     Acta Acustica, 2021.
%     
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/models/mckenzie2021.php

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

%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Stats
%   #Author: Thomas McKenzie


%%  DEFAULT OPTIONAL INPUTS

if nargin<8; lim = 0.05;   end
if nargin<7; plotFlag = 0; end
if nargin<6; w = 1;        end
if nargin<5; norm = [];    end

%%  FUNCTION SETTINGS

SPL = 75; % Sound Pressure Level
initInc = 0.2; % Initial Increment
scale = 0.4; 

% Factors control the impact of each perceptual modification made to the
% input signals. Selecting 0 for all factors will compute the same as an
% Actual Spectral Difference (A-B). Selecting 1 for each factor will
% compute the standard PSD. Factions are also accepted, which linearly vary
% the results inbetween an Actual and Perceptual Spectral Difference
% calculation. The factors should remain between -1 and 1 for 'meaningful'
% results.
pFac = 1; % Phons influence Factor
sFac = 1; % Sones influence Factor
erbFac = 1; % ERB influence Factor

% The Sones scale considers an increase in 10 Phons to be a doubling of
% Sones such that 40 Phons=1Sone. This factor allows that scale to be
% skewed within the PSD model such that an increase of sSkw phones results
% in a doubling of skewed Sones. 
sSkw = 10; % Sones skew factor

%%  CONVERT INPUTS TO FREQUENCY (dB) DOMAIN IF REQUIRED
    
switch domFlag
    % Time domain -> Frequency (dB) domain
    case 0 
        minfftSample = round(f.minFreq / f.fs * f.nfft)+1;
        maxfftSample = round(f.maxFreq / f.fs * f.nfft)+1;

        freq = (0 : f.fs/f.nfft : f.fs-(f.fs/f.nfft))';
        freq = freq(minfftSample:maxfftSample);

        A = 20*log10(abs(fft( A, f.nfft )));
        A = A(minfftSample:maxfftSample, :, :);

        B = 20*log10(abs(fft( B, f.nfft )));
        B = B(minfftSample:maxfftSample, :, :);

    % Frequency domain -> Frequency (dB) domain
    case 1
        A = 20*log10(abs(A));
        B = 20*log10(abs(B));
        freq = f(:);

    % Frequency (dB) domain
    case 2
        freq = f(:);

    otherwise
        error('Unrecognised integer input: domFlag. Must be 0, 1 or 2.');
end

%%  HANDLE INPUTS

sizeA = size(A);

% Map sample point weightings
if size(w) == 1
    w = repmat(w, [1, sizeA(2:end)]);
end

% Define maximum iterations / resolution limit
if lim >= 1
    maxIt = limit; % Maximum Iterations
    minRes = 0; % Minimum Resolution
elseif lim > 0
    maxIt = 100000;
    minRes = lim;
else
    error('Positive numeric value expected for input variable "limit"');    
end

offset = NaN(maxIt, 1);
err = NaN(maxIt, 1); % Error

%%  FUNCTION SETUP
    
% The code in PSD_Setup() is separated for convenience and computational
% efficiency. If PSD() is called repeatedly with the same frequency point
% analysis then is may be benificial to extract PSD_Setup(), store the
% output variables and adapt PSD() to accept the data as input variables.
[EL, ERBWeights] = setupPSD(sizeA, freq);
    
%%  FIND INPUT NORMALISATION
    
% Iterate to find perceptualy optimum value
if isempty(norm)

    % Initial comparason to find jump value
    [~, temp, jump] = calculatePSD(A, B, 0);
    temp_weighted = temp .* w;
    err(1) = sum(temp_weighted(:)) ./ sum(w(:));

    if plotFlag
        update = figure();title('Dataset normalisation'); hold on;
        plot(0, err(1), 'o'); hold on;
        drawnow
    end

    % For each iteration...
    for i = 2:maxIt

        % Find offset for this iteration
        switch i
            % Second comparason (No offset)
            case 2 
                offset(i) = jump;

            % Third comparason (Try out positive offset)
            case 3 
                offset(i) = jump + initInc;

            % All other comparasons
            otherwise 
                % If there was a reduction in perceptual error with respect
                % to the previous offset, carry on increasing the offset in
                % that direction (up or down)
                if err(i-1) < err(i-2) 
                    offset(i) = offset(i-1) + initInc;

                % If there was no improvement / an increase in error
                else 
                    % If this is the third iteration, switch the offset
                    % direction, but keep the offset resolution the same
                    if i == 4
                        initInc = initInc * -1;
                        offset(i) = jump + initInc;

                    % Otherwise, switch the offset direction and reduce the
                    % offset resolution
                    else
                        % If the resolution has got too fine, meaning
                        % reducing the offset resolution would make minimal
                        % difference to the PSD, quit
                        if abs(err(i-1)-err(i-2)) < minRes ||...
                           abs(err(i-2)-err(i-3)) < minRes
                            break;
                        else
                            initInc = initInc * -scale;
                            offset(i) = offset(i-1) + initInc;
                        end
                    end
                end
        end

        % Calculate PSD with calculated offset and save result
        [~, temp, ~] = calculatePSD(A, B, offset(i));
        temp_weighted = temp .* w;
        err(i) = sum(temp_weighted(:)) ./ sum(w(:));

        if plotFlag
            figure(update);
            plot(offset(i), err(i), 'o'); hold on;
            drawnow
        end

    end

    % Find optimal offset
    [~, minIndex] = min(err);
    norm = offset(minIndex);

    % Print result
    % fprintf('Optimum Normalisation Found: %d\n', offset_F);
end
    
%%  CALCULATE PERCEPTUAL SPECTRAL DIFFERENCE

[diff, avgDiff_ERB, ~] = calculatePSD(A, B, norm);
    
%%  INTERNAL FUNCTIONS
function [EL, ERBWeights] = setupPSD(sizeA, freq)

    % ISO 226 Declarations
    iso226SPL = zeros(30, 91);
    iso226Freq = zeros(30, 91);
    Y = zeros(length(freq), 91);
    EL = zeros(length(freq), 91);

    % For all ISO standardised listening levels (ISO 226: 0-90dB SPL)
    for l = 1:91
        % Save equal loudness contour
        [iso226SPL(1:29, l), iso226Freq(1:29, l)] = iso226(l-1);
        iso226Freq(30, l) = 20000;
        iso226SPL(30, l) = iso226SPL(1, l);

        % Fit curve to equal loudness contour
        iso226Fit = fit(iso226Freq(:, l), iso226SPL(:, l), 'pchip');

        % Interpolate to input frequency bins and remove equivalent 1KHz
        % loudness offset
        Y(:, l) = iso226Fit(freq) - (l - 1);

        % Save the offset required in dB to equate the loudness of any
        % frequency bin to that of 1KHz for a given absolute loudness
        % (0-90dB SPL)
        % ... A.K.A. flip the equal loudness contour 1KHz offset!
        EL(:, l) = -Y(:, l) ;
    end
    
    % Calculate ERB bandwidths
    ERB = 0.108.*freq + 24.7;

    % Calculate ERB weights and repeat for input matrix multiplication
    ERBWeights_temp = (1./ERB) ./ max(1./ERB);
    ERBWeights = repmat(ERBWeights_temp, [1, sizeA(2:end)]);
    
end

function [pDiff, avgPDiff_ERB, jump] = calculatePSD(A, B, norm)
%%  APPLY OFFSET AND NORMALISE INPUTS (75dB)

    B = B + norm;
    meanValue = mean([A(:); B(:)]);
    A = A + (SPL-meanValue);
    B = B + (SPL-meanValue);

%%  EQUAL LOUDNESS (PHONS)

    % Save matrices that select the correct loudness contour offsets to use
    % for each frequency bin for each input spectrum. This is calculated by
    % linear interpolation of the integer rounded input matricies with a
    % maximum value of 90 and a minimum value of 0.
    
    % Round to hi/lo intergers and save fractional component
    LC_A_hi = min(ceil(A),  90); LC_A_hi = max(LC_A_hi, 0);
    LC_A_lo = min(floor(A), 90); LC_A_lo = max(LC_A_lo, 0);
    A_frac = rem(A, 1);
    LC_B_hi = min(ceil(B),  90); LC_B_hi = max(LC_B_hi, 0);
    LC_B_lo = min(floor(B),  90); LC_B_lo = max(LC_B_lo, 0);
    B_frac = rem(B, 1);

    % find loudness contour offsets for integer values
    EL_A_hi = EL(sub2ind(size(EL),...
                         repmat((1:sizeA(1))', [1, sizeA(2:end)]),...
                         LC_A_hi+1));
    EL_A_lo = EL(sub2ind(size(EL),...
                         repmat((1:sizeA(1))', [1, sizeA(2:end)]),...
                         LC_A_lo+1));
    EL_B_hi = EL(sub2ind(size(EL),...
                         repmat((1:sizeA(1))', [1, sizeA(2:end)]),...
                         LC_B_hi+1));
    EL_B_lo = EL(sub2ind(size(EL),...
                         repmat((1:sizeA(1))', [1, sizeA(2:end)]),...
                         LC_B_lo+1));

    % Interpoloate between integer based loudness contour offsets
    EL_A = EL_A_lo + (A_frac.*(EL_A_hi - EL_A_lo));
    EL_B = EL_B_lo + (B_frac.*(EL_B_hi - EL_B_lo));
    
    % Account for equal loudness / convert input data to phones scale
    A_EL = A + EL_A;
    B_EL = B + EL_B;
    
    % Apply Factor
    A_EL = ((1-pFac) .* A) + (pFac .* A_EL);
    B_EL = ((1-pFac) .* B) + (pFac .* B_EL); 

%%  CONVERT TO SONES

    A_EL_Sones = 2.^((A_EL-40)/sSkw);
    B_EL_Sones = 2.^((B_EL-40)/sSkw);
    
    % Apply Factor
    A_EL_Sones = ((1-sFac) .* A_EL) + (sFac .* A_EL_Sones);
    B_EL_Sones = ((1-sFac) .* B_EL) + (sFac .* B_EL_Sones);

%%  FIND PERCEPTUAL AVERAGE DIFFERENCE WITH ERB WEIGHTING

    pDiff = (B_EL_Sones - A_EL_Sones);
    avgPDiff_ERB = sum(ERBWeights .* abs(pDiff)) ./ sum(ERBWeights);

    % Apply Factor
    avgPDiff_ERB = ((1-erbFac) .* mean(abs(pDiff))) + (erbFac .* avgPDiff_ERB);
    
    % Jump provides a rough idea of the level difference required to make
    % to two inputs sound more similar
    A_level = sum(ERBWeights .* (A_EL_Sones)) ./ sum(ERBWeights);
    B_level = sum(ERBWeights .* (B_EL_Sones)) ./ sum(ERBWeights);
    jump = mean(A_level(:)) - mean(B_level(:));
        
end

end

function [spl, freq] = iso226(phon)
% Equal Loudness Curves from ISO 226

f = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 ...
     1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500];

af = [0.532 0.506 0.480 0.455 0.432 0.409 0.387 0.367 0.349 0.330 0.315 ...
      0.301 0.288 0.276 0.267 0.259 0.253 0.250 0.246 0.244 0.243 0.243 ...
      0.243 0.242 0.242 0.245 0.254 0.271 0.301];

Lu = [-31.6 -27.2 -23.0 -19.1 -15.9 -13.0 -10.3 -8.1 -6.2 -4.5 -3.1 ...
       -2.0  -1.1  -0.4   0.0   0.3   0.5   0.0 -2.7 -4.1 -1.0  1.7 ...
        2.5   1.2  -2.1  -7.1 -11.2 -10.7  -3.1];

Tf = [ 78.5  68.7  59.5  51.1  44.0  37.5  31.5  26.5  22.1  17.9  14.4 ...
       11.4   8.6   6.2   4.4   3.0   2.2   2.4   3.5   1.7  -1.3  -4.2 ...
       -6.0  -5.4  -1.5   6.0  12.6  13.9  12.3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%Error Trapping
if((phon < 0) | (phon > 90))
    disp('Phon value out of bounds!')
    spl = 0;
    freq = 0;
else
    %Setup user-defined values for equation
    Ln = phon;

    %Deriving sound pressure level from loudness level (iso226 sect 4.1)
    Af=4.47E-3 * (10.^(0.025*Ln) - 1.15) + (0.4*10.^(((Tf+Lu)/10)-9 )).^af;
    Lp=((10./af).*log10(Af)) - Lu + 94;

    %Return user data
    spl = Lp;  
    freq = f;
end
end

