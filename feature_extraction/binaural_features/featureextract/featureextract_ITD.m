function [ ITDs, xc, rHat ] = featureextract_ITD(FEXT, in)
%FEATUREEXTRACT_ITD get ITDs and optionally crosscorrelations

% input is expected to be (blocksize x N x F x M)
if ndims(in)<4
    error('Input wrong shape.');
end

nLags = length(FEXT.lagsITD);

N = size(in, 2);
F = size(in, 3);

rHat = zeros(N, F);
ITDs = zeros(N, F);
xc = zeros(nLags, N, F);
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
    ITDs(:, ii) = (FEXT.lagsITD(bM) + delta.')/FEXT.fs;
    xc(:, :, ii) = xcorr;
end

end

