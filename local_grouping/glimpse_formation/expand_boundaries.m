function [expanded_boundaries] = expand_boundaries(boundary_map)
%BOUNDARY_EXPANSION expand boundaries of an already existing binary matrix
%or image by adding boundary pixels in every direction
%   
% go through every pixel of the map
% if any direct neighbor is a 1 , set the current pixel to 1

expanded_boundaries = zeros(size(boundary_map));

for aa = 2:(size(boundary_map,1)-1)
    for bb = 2:(size(boundary_map,2)-1)
        testsum = boundary_map(aa-1,bb)+boundary_map(aa+1,bb)+boundary_map(aa,bb-1)+boundary_map(aa,bb+1)+boundary_map(aa,bb);
        if testsum >= 1
            expanded_boundaries(aa,bb)=1;
        else 
            expanded_boundaries(aa,bb)=0;
        end
    end
end

for aa = 2:(size(boundary_map,1)-1)
    testsum1    = boundary_map(aa-1,1)+boundary_map(aa+1,1)+boundary_map(aa,2)+boundary_map(aa,1);
    if testsum1 >= 1
        expanded_boundaries(aa,1) = 1;
    end
    
    testsumEnd  = boundary_map(aa-1,end)+boundary_map(aa+1,end)+boundary_map(aa,(end-1))+boundary_map(aa,end);
    if testsumEnd >= 1
        expanded_boundaries(aa,end) = 1;
    end
end



end

