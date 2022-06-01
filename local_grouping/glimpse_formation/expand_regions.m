function [expanded_regions] = expand_regions(regions_map, num_regions, mode, boundary_diffs, current_bound)
%EXPAND_REGIONS expand regions in an already existing binary matrix
%or image by adding pixels in every direction of the existing region
%   
% INPUTS: 
%   regions_map     - binary matrix/image, regions are labeled with ones
%   num_regions     - number of regions (clusters of ones) in regions_map
%   mode            - can be set to 'undirected' or 'directed'
%                     'undirected' means, go through every pixel of the map
%                     if any direct neighbor is a 1 , set the current pixel to 1
%                     to use 'directed', the input boundary_diffs is needed
%                     which indicates different stages of boundary
%                     expansion. In this mode, the regions are expanded in
%                     a way that they first fill those pixels that were
%                     last occupied by a boundary expansion
%   boundary_diffs  - indicates different stages of boundary expansion
%                     (shows the iteratively built up boundaries)
%   current_bound   - the iteration of boundary expansion that's currently
%                     considered (boundary that regions are allowed to grow
%                     into)
% OUTPUTS:
%   expanded_regions- same dimensions as regions_map but with expanded
%                     regions


new_regions_map = zeros(size(regions_map));

if strcmp(mode,'undirected') == 1  % the regions grow one pixel in every direction, starting from inital regions map
    
    for aa = 1:num_regions
        
        [rows, cols] = find(regions_map==aa);
        
        % expand every region by one pixel in each direction
        
        next1 = [rows(:)-1, cols(:)];
        next2 = [rows(:)+1, cols(:)];
        next3 = [rows(:), cols(:)-1];
        next4 = [rows(:), cols(:)+1];
        
        next5 = [rows(:)-1, cols(:)-1];
        next6 = [rows(:)+1, cols(:)+1];
        next7 = [rows(:)-1, cols(:)+1];
        next8 = [rows(:)+1, cols(:)-1];
        
        all_neighbors = [next1; next2; next3; next4; next5; next6; next7; next8];
        unique_inds = unique(all_neighbors,'rows');
        
        [rows_rem1, ~] = find(unique_inds == 0);
        unique_inds(rows_rem1,:) = [];    % remove all rows that contain a zero
        [rows_rem2, ~] = find(unique_inds(:,1) > size(regions_map,1));
        unique_inds(rows_rem2,:) = [];    % remove all rows that contain indices which exceed the matrix dimensions of the regions_map
        [rows_rem3, ~] = find(unique_inds(:,2) > size(regions_map,2));
        unique_inds(rows_rem3,:) = [];

        % exclude all rows that contain a zero, a negative number, or an index out of bounds of the regions_map
        
        for bb = 1:size(unique_inds,1)
            new_regions_map(unique_inds(bb,1),unique_inds(bb,2))=aa;
        end
        
    end
elseif strcmp(mode,'directed') == 1
    % needs to grow into boundary_diffs
    
    for aa = 1:num_regions
        
        [rows, cols] = find(regions_map==aa);
        
        % expand every region by one pixel in each direction
        
        next1 = [rows(:)-1, cols(:)];
        next2 = [rows(:)+1, cols(:)];
        next3 = [rows(:), cols(:)-1];
        next4 = [rows(:), cols(:)+1];
        
        next5 = [rows(:)-1, cols(:)-1];
        next6 = [rows(:)+1, cols(:)+1];
        next7 = [rows(:)-1, cols(:)+1];
        next8 = [rows(:)+1, cols(:)-1];
        
        all_neighbors = [next1; next2; next3; next4; next5; next6; next7; next8];
        unique_inds = unique(all_neighbors,'rows');
        
        [rows_rem1, ~] = find(unique_inds == 0);
        unique_inds(rows_rem1,:) = [];    % remove all rows that contain a zero
        [rows_rem2, ~] = find(unique_inds(:,1) > size(regions_map,1));
        unique_inds(rows_rem2,:) = [];    % remove all rows that contain indices which exceed the matrix dimensions of the regions_map
        [rows_rem3, ~] = find(unique_inds(:,2) > size(regions_map,2));
        unique_inds(rows_rem3,:) = [];

        % exclude all rows that contain a zero, a negative number, or an index out of bounds of the regions_map
        
        for bb = 1:size(unique_inds,1)
%             if boundary_diffs(unique_inds(bb,1),unique_inds(bb,2))>=current_bound  || regions_map(unique_inds(bb,1),unique_inds(bb,2))==aa
            %if 0 >= boundary_diffs(unique_inds(bb,1),unique_inds(bb,2))>=current_bound  || regions_map(unique_inds(bb,1),unique_inds(bb,2))==aa
            if boundary_diffs(unique_inds(bb,1),unique_inds(bb,2))==current_bound  || regions_map(unique_inds(bb,1),unique_inds(bb,2))==aa
                % replace the pixel 
                new_regions_map(unique_inds(bb,1),unique_inds(bb,2))=aa;
                
            elseif boundary_diffs(unique_inds(bb,1),unique_inds(bb,2))~=current_bound && regions_map(unique_inds(bb,1),unique_inds(bb,2))~=aa
                % do not replace the pixel but leave it a zero
                new_regions_map(unique_inds(bb,1),unique_inds(bb,2))=0;
            end
            
        end
    end
    
else
    error('define mode for region expansion')
end

expanded_regions = new_regions_map;
% fix, what happens, when 2 regions start overlapping each other
% fix what happens, when some pixels remain empty


% check if this belongs to the next ring


end

