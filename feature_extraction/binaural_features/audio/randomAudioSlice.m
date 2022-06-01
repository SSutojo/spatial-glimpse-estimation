function aObj = randomAudioSlice(aObj, nSamples)
%randomAudioSlice grab random slice of nSamples of input audio

%   ***********************************************************************

offset = randi(aObj.nSamples-nSamples-1);
aObj.data = aObj.data((1:nSamples)+offset, :);
if isfield(aObj, 'label')
    aObj.label = [aObj.label ' (slice)'];
end
aObj = updateAudio(aObj);

%   ***********************************************************************