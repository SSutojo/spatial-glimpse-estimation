function [ noise ] = spat_genDiffNoisePINK(SPAT, L)
%SPAT_GENDIFFNOISE generate spatially diffuse pink noise of length L

nDirs = size(SPAT.diffuseDirs, 2);
noise = zeros(L, SPAT.nChans);

beta = -1;
freqs = [linspace(0,.5,floor(L/2)),fliplr(linspace(0,.5,ceil(L/2)))]';
freqs = (freqs.^2).^(beta/2);
freqs(freqs==inf) = 0;
phase_vals = rand(L, 1);


for n=1:nDirs
%     noiseSSN = fftfilt(b, randn(L, 1));
    noiseSSN = real(ifft(sqrt(freqs).*(cos(2*pi*phase_vals)+1i*sin(2*pi*phase_vals))));
    noise = noise + spatialize_bin( SPAT, noiseSSN, ...
        SPAT.diffuseDirs(1, n), SPAT.diffuseDirs(2, n));
end

noise = genAudio(noise, SPAT.fs, 'Diffuse Noise');

end


% psd überprüfen mit [Pxx,W] = pwelch(noiseSSN,512,256,512);
% und semilogx(W,10*log10(Pxx))