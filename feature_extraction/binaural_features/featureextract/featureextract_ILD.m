function [ ILDs ] = featureextract_ILD(FEXT, in)
%FEATUREEXTRACT_ILD Extract ILD from input (chans 1&2 only)

% input is expected to be (blocksize x N x F x M)
if ndims(in)<4
    error('Input wrong shape.');
end

bmEnergy = shiftdim(sum(power(in(:,:,:,1:2),2))/FEXT.blocksize);
ILDs = 10 * (log10(bmEnergy(:,:,2) + eps) -  log10(bmEnergy(:,:,1) + eps));
if isequal(FEXT.HaircellModel,'roman')
    % Reverse effect of sqrt()
    ILDs = 2 * ILDs;
end

end

