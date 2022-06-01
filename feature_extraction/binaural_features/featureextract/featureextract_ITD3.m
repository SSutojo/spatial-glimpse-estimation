function [ ITDs, rHat ] = featureextract_ITD3(FEXT, in)
%FEATUREEXTRACT_ITD3 get ITDs of a 4-channel signal

% input is expected to be (blocksize x N x F x M)
if ndims(in)<4
    error('Input wrong shape.');
end

nLags = length(FEXT.lagsITD);

N = size(in, 2);
F = size(in, 3);

maxDelayx = ceil(FEXT.maxTDOAsec * FEXT.fs)+1;
lagsx = (-maxDelayx:1:maxDelayx);
nLagsx = length(lagsx);

rHat = zeros(N, F);
ITDs = zeros(N, F, 2);
I1 = zeros(N, F);
I2 = zeros(N, F);
%I3 = zeros(N, F);
for ii=1:F
    xcorr = xcorrNorm(in(:,:, ii, 1), in(:,:, ii, 2), ...
        FEXT.maxDelay, true, true);
    % Maximum search
    [~,b,~,~] = findLocalPeaks(xcorr,2);
    % Warp global maximum indices
    bM = mod(b-1,nLags)+1;
    % Parabolic interpolation
    [delta,rHat(:, ii)] = interpolateParabolic(xcorr,bM);
    % Calculate ITD
    I1(:, ii) = (FEXT.lagsITD(bM) + delta.')/FEXT.fs;

    % ITD left side mics
    xcorr = xcorrNorm(in(:,:, ii, 1), in(:,:, ii, 3), ...
        maxDelayx, true, true);
    [~,b,~,~] = findLocalPeaks(xcorr,2);
    bM = mod(b-1,nLagsx)+1;
    [delta,rHat_l] = interpolateParabolic(xcorr,bM);
    %I2(:, ii) = (lagsx(bM) + delta.')/FEXT.fs;
    lags_l = (lagsx(bM) + delta.')/FEXT.fs;
    rHat_l((bM==1)|(bM==nLags+1)) = 0;
    
    % ITD right side mics
    xcorr = xcorrNorm(in(:,:, ii, 2), in(:,:, ii, 4), ...
        maxDelayx, true, true);
    [~,b,~,~] = findLocalPeaks(xcorr,2);
    bM = mod(b-1,nLagsx)+1;
    [delta,rHat_r] = interpolateParabolic(xcorr,bM);
    %I3(:, ii) = (lagsx(bM) + delta.')/FEXT.fs;
    lags_r = (lagsx(bM) + delta.')/FEXT.fs;
    rHat_r((bM==1)|(bM==nLags+1)) = 0;
    
    l_or_r = rHat_l > rHat_r;
    I2(l_or_r, ii) = lags_l(l_or_r);
    I2(~l_or_r, ii) = lags_r(~l_or_r);
end

ITDs(:, :, 1) = I1;
ITDs(:, :, 2) = I2;
%ITDs(:, :, 3) = I3;

end

