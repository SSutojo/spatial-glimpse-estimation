function [ noise ] = spat_genDiffNoise(SPAT, L, b)
%SPAT_GENDIFFNOISE generate diffuse noise of length L, shaped by filter b

nDirs = size(SPAT.diffuseDirs, 2);
noise = zeros(L, SPAT.nChans);

for n=1:nDirs
    noiseSSN = fftfilt(b, randn(L, 1));
    noise = noise + spatialize( SPAT, noiseSSN, ...
        SPAT.diffuseDirs(1, n), SPAT.diffuseDirs(2, n));
end

noise = genAudio(noise, SPAT.fs, 'Diffuse Noise');

end