function [ ILDs, ITDs, ICs ] = featureextract_FFTbased( FEXT, insig )
%FEATUREEXTRACT_FFTBASED calculate FFT based features

insize = size(insig.data);
L = insize(1);
M = insize(2);

numblocks = floor(L./FEXT.hopsize)-1;
insig_blocked = zeros(FEXT.blocksize, numblocks, M);

for m=1:M
    insig_blocked(:,:,m) = frameData(insig.data(:,m), ...
        FEXT.blocksize, FEXT.hopsize, hann(FEXT.blocksize));
end

STFTsig = fft(insig_blocked, FEXT.fftsize, 1); % (FsizexNxM)
ILDs = zeros(numblocks, FEXT.nFilts);
ITDs = zeros(numblocks, FEXT.nFilts);
ICs = zeros(numblocks, FEXT.nFilts);

% calculate ILD, IPD, IC (Breebart 2005)
for f=1:FEXT.nFilts
    fbin = FEXT.fb_ss(1,f);
    lbin = FEXT.fb_ss(2,f);
    X1 = STFTsig(fbin:lbin, :, 2);
    X2 = STFTsig(fbin:lbin, :, 1);
    X1ms = sum(abs(X1).^2, 1);
    X2ms = sum(abs(X2).^2, 1);
    cross = sum(X2.*conj(X1), 1);
    ILDs(:,f) = 10*log10(X1ms ./ X2ms);
    ITDs(:,f) = angle(cross);
    ICs(:, f) = abs(cross) ./ sqrt(X1ms .* X2ms);
end

end

