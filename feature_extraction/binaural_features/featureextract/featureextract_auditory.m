function [ afout, FEXT ] = featureextract_auditory(FEXT, x, bNeural)
%FEATUREEXTRACT_AUDITORY Apply auditory frontend (FULL SIG ONLY)

narginchk(2, 3);
if ~isAudio(x)
    error('Input is not Audio object.');
end

if nargin<3
    bNeural = true;
end

x = x.data;
L = size(x, 1);
afout = zeros(L, FEXT.nFilts, FEXT.nChans);

for ii = 1:FEXT.nChans
    [afout(:, :, ii), FEXT] = FEXT.gtfb(x(:,ii), ii, FEXT);
end

md = max(FEXT.delays);
for mm = 1:FEXT.nFilts
    d = md-FEXT.delays(mm);
    if d>0
        temp = afout(L-d+1:L, mm, :);
        afout(:, mm, :) = [FEXT.delaymem{mm}; afout( 1:L-d, mm, :)];
        FEXT.delaymem{mm} = temp;
    end
end

if bNeural
    % apply neural transduction model
    afout = FEXT.haircell(afout); %sqrt(max(afout, 0));
end

end
