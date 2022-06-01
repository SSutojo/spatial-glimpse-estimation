function [ outsig ] = featureextract_blockFFT( FEXT, insig )
%FEATUREEXTRACT_BLOCK split input into short-time blocks

insize = size(insig.data);
L = insize(1);
F = FEXT.fftsize/2+1;
M = insize(2);

numblocks = floor(L./FEXT.hopsize)-1;
insig_blocked = zeros(FEXT.blocksize, numblocks, M);

for m=1:M
    insig_blocked(:,:,m) = frameData(insig.data(:,m), ...
        FEXT.blocksize, FEXT.hopsize, hann(FEXT.blocksize));
end

outsig = fft(insig_blocked, FEXT.fftsize, 1);
outsig = permute(outsig(1:F, :, :), [4 2 1 3]);

end

