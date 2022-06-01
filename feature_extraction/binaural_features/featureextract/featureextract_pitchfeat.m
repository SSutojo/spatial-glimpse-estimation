function [ fout ] = featureextract_pitchfeat( FEXT, x )
%FEATUREEXTRACT_PITCHFEAT Extract Pitch features

narginchk(2, 3);
if ~isAudio(x)
    error('Input is not Audio object.');
end

x = x.data;
L = size(x, 1);
afout = zeros(L, FEXT.nFilts, FEXT.nChans);
envout = zeros(L, FEXT.nFilts, FEXT.nChans);

% gammatone filtering
for ii = 1:FEXT.nChans
    [afout(:, :, ii), envout(:, :, ii)] = ...
        gammatoneMEX(x(:,ii), FEXT.fs, FEXT.cfs(1), FEXT.cfs(end), ...
        [FEXT.nFilts 1:FEXT.nFilts], 0, 1, 0);
    % delay to align
    afout(1:L-max(FEXT.delays),ii,:) = ...
        afout(FEXT.delays(ii)+(1:L-max(FEXT.delays)),ii,:);
    envout(1:L-max(FEXT.delays),ii,:) = ...
        envout(FEXT.delays(ii)+(1:L-max(FEXT.delays)),ii,:);
end

% if cf>1.5kHz, use envelopes
afout(:, FEXT.cfs>1500, :) = envout(:, FEXT.cfs>1500, :);

% blocking
blafdata = featureextract_block(FEXT, afout);

slag = 38;
elag = 228;
% calculate NAC and CFR - on both channels
for j = 1:size(blafdata, 2)          % iterate over frames
    for k = 1:size(blafdata, 3)      % iterate over auditory channels
        for m = 1:size(blafdata, 4)  % iterate over microphone channels
            nacdata(:, j, k, m) = NAC(blafdata(:, j, k, m), slag, elag);
            cfrdata(:, j, k, m) = CFR(blafdata(:, j, k, m), slag, elag);
        end
    end
end
PD = max(0.01, nacdata.*cfrdata);
PD = mean(PD, 4);
PD = mean(PD, 3);

end

function out = NAC(in, s, e)

L = length(in)-e+1;

x1 = in(1:L);
ex1 = sum(x1.^2);
out = zeros(1, e-s+1);
for p=s:e
    x2 = in(p:p+L-1);
    ex2 = sum(x2.^2);
    out(p-s+1) = sum(x1.*x2)/(ex1*ex2);
end

end

function out = CFR(in, s, e)

L = length(in)-e+1;

x1 = in(1:L);
out = zeros(1, e-s+1);
for p=s:e
    x2 = in(p:p+L-1);
    sumx1x2 = x1+x2;
    diffx1x2 = x1-x2;
    out(p-s+1) = sum(sumx1x2.^2)/sum(diffx1x2.^2);
end

end
