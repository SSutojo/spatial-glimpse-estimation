function [ lSNR ] = featureextract_localSNR(FEXT, speech, noise)
%LOCALSNR find t/f local SNR

if FEXT.fb~='FFT'
    s = featureextract_auditory(FEXT, speech, false);
    n = featureextract_auditory(FEXT, noise, false);
    
    sb = featureextract_block( FEXT, s );
    nb = featureextract_block( FEXT, n );
    
    sb = shiftdim(sum(abs(sb).^2, 1)/FEXT.blocksize, 1);
    nb = shiftdim(sum(abs(nb).^2, 1)/FEXT.blocksize, 1);
    
    sb = sum(sb, 3);
    nb = sum(nb, 3);
    
    lSNR = 10*log10(sb./nb);
else
    sb = featureextract_block( FEXT, speech.data );
    nb = featureextract_block( FEXT, noise.data );

    sbSTFT = sum(fft(sb, FEXT.fftsize, 1), 3); % (FsizexN)
    nbSTFT = sum(fft(nb, FEXT.fftsize, 1), 3);
    
    lSNR = zeros(size(sbSTFT, 2), FEXT.nFilts);
    for f=1:FEXT.nFilts
        fbin = FEXT.fb_ss(1,f);
        lbin = FEXT.fb_ss(2,f);
        s = sbSTFT(fbin:lbin, :);
        n = nbSTFT(fbin:lbin, :);
        s = sum(abs(s).^2, 1);
        n = sum(abs(n).^2, 1);
        lSNR(:, f) = 10*log10(s./n);
    end

end
    

end

