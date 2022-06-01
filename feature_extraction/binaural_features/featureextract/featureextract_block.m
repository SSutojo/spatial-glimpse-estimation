function [ outsig ] = featureextract_block( FEXT, insig )
%FEATUREEXTRACT_BLOCK split input into short-time blocks

insize = size(insig);
L = insize(1);
F = insize(2);
if length(insize)>2
    M = insize(3);
else
    M = 1;
end

% numblocks = floor(L./FEXT.hopsize)-1;
numblocks = floor((L-(FEXT.blocksize-FEXT.hopsize))/FEXT.hopsize);
outsig = zeros(FEXT.blocksize, numblocks, F, M);

for f=1:F
    for m=1:M
        outsig(:,:,f,m) = frameData(insig(:,f,m), ...
            FEXT.blocksize, FEXT.hopsize, FEXT.blkwindow);
    end
end

end

