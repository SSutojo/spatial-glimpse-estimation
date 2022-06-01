function [G, num_regions, seed_image, TI] = regiongrow(F, S, T)
%REGIONGROW Perform segmentation by region growing
% Function uses an image and forms segments of similar colour (or intensity) e.g. "grows regions" in this image
% (taken from Gonzales, Woods, Eddins "Digital Image Processing
% Using MATLAB", chapter 10, page 423)
% The regions start to grow from "seeds" that can be manually placed in the
% image or 
%
% INPUTS:   F:  is an image to be segmented. array which represents the map in which the regions are
%               supposed to be grown (e.g. within contours in that map)
%           S:  can be an array(the same size as F) with a 1 at the coordinates of
%               every seed point and 0s elsewhere. S can also be a single seed value.
%               If S is an array the same size as F, it contains a 1 at
%               each seed point.
%               If S is a scalar it defines an intensity value such that all the points
%               in F with that intensity value become seed points
%           T:  can be an array (the same size as F) containing a threshold 
%               value for each pixel in F. T can also be a scalar, in which case it
%               becomes a global threshold.
%
% OUTPUTS:  G:  is the result of region growing, with each region
%               labeled by a different intergerm 
%           num_regions: is the number of regions, 
%           seed_image: is the final seed image used by the algorithm 
%           TI: is the image consisting of the pixels in F that satisfied the threshold test.



F = double(F);
% if S is a scalar, obtain the seed image.
if numel(S) == 1
    seed_image = F == S;
    S1 = S;
else 
    %S is an array. Eliminate duplicate, connected seed locations to reduce
    %the number of loop executions in the following sections of code
    seed_image = bwmorph(S, 'shrink', Inf);
    J = find(seed_image);
    S1 = F(J); % array of seed values
end

TI = false(size(F));
for k = 1:length(S1)
    seedvalue = S1(k);
    S = abs(F - seedvalue) <= T;
    TI = TI | S;
end

% use function imreconstruct with seed_image as the marker image to onbtain
% the regions corresponding to each seed in S. Function bwlabel assigns a
% different integer to each connected region. 
[G, num_regions] = bwlabel(imreconstruct(seed_image, TI));


end


