function [regions_est,all_region_maps,all_boundary_maps,all_merged_boundary_maps] = SoftBoundaryExp(contrast_map,thresh_vals,Nexpansions)
%SOFTBOUNDARYREGIONSEXP repeat the Boundary and Region Expansion
%with shifting contrast thresholds 
%
% IN:
%       contrast_map - matrix with contrast values (0-1) for each transition between T-F units
%                      Dims: (number_time_frames*2-1) X (number_subbands*2-1)
%                      placeholders for T-F units are zeros
%       thresh_vals  - vector with different contrast thresholds to apply to the contrast_map and get different boundar_maps 
%                      (the lower this threshold, the more conservative)
%                      start with the lowest threshold to get the most
%                      conservative boundary map
%       Nexpansions  - number of time that the initial boundaries are expanded
% 
% OUT:   
%       regions_est  - map the indicates estimated regions (as connected
%                      regions labeled with the same integer) and
%                      boundaries between these regions (labeled with
%                      zeros)
%                      Dims: (number_time_frames*2-1) X (number_subbands*2-1)
%       all_region_maps - regions_est from each iteration
%                      Dims: (number_time_frames*2-1) X (number_subbands*2-1) X numel(thresh_vals)
%       all_boundary map - all initial boundary maps for each iteration, directly after applying the
%                      current contrast threshold to the contrast map
%                      Dims: (number_time_frames*2-1) X (number_subbands*2-1) X numel(thresh_vals)
%       all_merged_boundary_maps - boundary maps (after boundary expansion
%                      and merging with previous maps) for each iteration
%                      can be used to visualize how the seed glimpses change in each iteration
%                      Dims: (number_time_frames*2-1) X (number_subbands*2-1) X numel(thresh_vals)


%%
previous_regions            = zeros(size(contrast_map,1),size(contrast_map,2));
boundary_IDs                = [1:Nexpansions].*(-1);

all_region_maps             = zeros(size(previous_regions,1),size(previous_regions,2),length(thresh_vals));
all_boundary_maps           = zeros(size(previous_regions,1),size(previous_regions,2),length(thresh_vals));
all_merged_boundary_maps    = zeros(size(previous_regions,1),size(previous_regions,2),length(thresh_vals));

% loop over different contrast thresholds
for aa = 1:length(thresh_vals)  % start with a low criterion and end with high criterion
    % get initial boundary map, by setting current threshold
    initial_boundary_map        = (contrast_map >= thresh_vals(aa));
    aa_boundary_maps            = zeros(size(initial_boundary_map,1),size(initial_boundary_map,2),(Nexpansions+1));
    aa_boundary_maps(:,:,1)     = initial_boundary_map;
    
    % expand the boundaries ('dilation' of boundaries)
    temp_boundary_map           = initial_boundary_map;
    for bb = 1:Nexpansions
        boundary_map_bb             = expand_boundaries(temp_boundary_map);
        aa_boundary_maps(:,:,bb+1)  = boundary_map_bb;
        temp_boundary_map           = boundary_map_bb;
    end
    largest_boundary_map        = aa_boundary_maps(:,:,end);
    
    % get boundary_diffs by stepwise subtraction of the smaller boundary maps
    aa_boundary_diffs              = largest_boundary_map.*(-1);  % subtract the smaller boundary maps
    for bb = 1:Nexpansions
        aa_boundary_diffs          = aa_boundary_diffs - aa_boundary_maps(:,:,(end-bb));
    end
    
    % merge the current boundary map and previously formed glimpses by occuppying zero regions in the boundary map that are already occupied by previously grown regions  
    aa_merged               = largest_boundary_map + previous_regions;
    aa_merged_boundaries    = (aa_merged >= 0.5); % leaves only regions open where no boundary and no previous region is, these are new seeds
    
    % use regiongrow in the remaining zeros, these are new seeds
    [added_seeds, ~, ~, ~] = regiongrow(aa_merged_boundaries,0,0.2);  % regiongrow in the remaining zeros
    
    
    % region expansion ('erosion' of boundaries)
    aa_region_maps          = zeros(size(temp_boundary_map,1),size(temp_boundary_map,2),Nexpansions+1);
    Nseeds                  = max(max(added_seeds));
    aa_region_maps(:,:,1)   = added_seeds;
    temp_region_map = added_seeds;

    for bb = 1:Nexpansions     
        curr_boundary           = boundary_IDs(bb);
        bb_regions_map          = expand_regions(temp_region_map, Nseeds,'directed', aa_boundary_diffs, curr_boundary);
        temp_region_map         = bb_regions_map;
        aa_region_maps(:,:,bb+1)= bb_regions_map;
    end
    
    % merge the newly grown regions with regions from previous iteration
    added_regions                   = aa_region_maps(:,:,end);
    AS                              = (added_regions >= 0.5);
    added_regions                   = (added_regions + max(max(previous_regions))).* AS;
    largest_region_map              = previous_regions + added_regions;
    
    all_region_maps(:,:,aa)         = largest_region_map; % keep the current glimpse map
    all_boundary_maps(:,:,aa)       = initial_boundary_map; 
    all_merged_boundary_maps(:,:,aa)= aa_merged_boundaries;
    previous_regions                = largest_region_map;  % save for next iteration   
    
end

regions_est = all_region_maps(:,:,end);


end

