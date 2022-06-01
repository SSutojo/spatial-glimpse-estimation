function [azmMat] = drawAzimuths(sepAngle,nSources, nMixes)
%DRAWAZIMUTHS generate azimuth values for each individual source in the mix
%for a given number of mixtures (nMixes)
% 
% INPUTS
%   sepAngle : desired separation angle between the sources
%   nSources : number of sources for which an azimuth is drawn 
%   nMixes   : number of mixes 
% 
% OUTPUT
%   azmMat   : matrix with azm per source and per mix. Dims: (nMixes X nSources)

%% check if the separation angle is too big to fit all sources
neededAngle = nSources * sepAngle;

if neededAngle>180
    error('separation angle and number of sources do not allow equal spacing in azimuth plane')
end

%% azimuth position of the first source
start_pos       = -90:5:(90-(sepAngle*(nSources-1))); % possible start positions
idx_vec         = randi(numel(start_pos),nMixes,1);

%% azimuth positions of other sources according to the desired separation angle
azmMat          = zeros(nMixes, nSources);
for aa = 1:nSources
    % draw an azimuth for each source
    azmMat(:,aa) = start_pos(idx_vec) + (aa-1)*sepAngle;
end



end

