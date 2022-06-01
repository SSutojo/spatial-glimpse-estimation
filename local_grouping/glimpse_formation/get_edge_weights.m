function [ edge_weights, cluster_size , cluster_element_positions ] = get_edge_weights( weights_and_labels_map , current_label )
%GET_EDGE_WEIGHTS get the edge weights within a cluster in which all
%t-f-units have the same label
%
% INPUTS:
%       weights_and_labels_map - contains the integers indicating the formed
%                   superpixels and the contrast values at the
%                   transitions. Dimensions: (time_frames*2-1 , frequency_bands*2-1)
%       current_label - the currently regarded glimpse (each TF unit that's
%                   labeled with this integer belongs to the same
%                   glimpse)
% 
% OUTPUTS:
%       edge_weights - vector that contains all weight (i.e. contrast
%                   values) that are within the glimpse
%       cluster_size - size of the glimpse in nr. of TF units
%       cluster_element_positions - time frame indices and filterband
%                   indices of the TF units belonging to this glimpse 


edge_weight_positions = [];


% find indices of the t-f-units with the current label in weights_and_labels_map

[rows, cols] = find(weights_and_labels_map==current_label);
cluster_size = numel(rows);   % number of t-f units with this label is the size of the cluster
cluster_elements = [rows, cols]; % coordinates of the cluster elements (x,y)
if cluster_size == 1
    edge_weights = 0;
    
elseif cluster_size < 1
    error('label does not exist')
elseif cluster_size > 1
    
    % find the t-f-units with the same x-coordinate (rows)
    U = unique(rows);
    repeated_x_coords = U(histc(rows,unique(rows))>1); % the x coordinates that appear more than once in rows
    for aa=1:length(repeated_x_coords) % go through every x coordinate that appears more than once
        x_inds = find(cluster_elements(:,1)==repeated_x_coords(aa));
        same_row_elements = cluster_elements(x_inds,:);        
        % sort them according to their y coordinate
        same_row_elements = sort(same_row_elements);
        % size(same_row_elements,1) = number of elements
        if size(same_row_elements,1)>1          
            for bb = 2:size(same_row_elements,1)
                dist = abs(same_row_elements(bb,:)-same_row_elements(bb-1,:)); % subtract the y-coordinates of neighboring elements
                distance = sum(abs(dist));
                if distance == 2 % if the difference is exactly 2 then they are neighbors
                    % if they are neighbors, the edge weight connecting them should be considered
                    curr_x_weight=(same_row_elements(bb,:)+same_row_elements(bb-1,:))/2; % coordinates of the edge_weight
                    edge_weight_positions = [edge_weight_positions ; curr_x_weight];
                end
            end
        end
    end
    % do the same for same column elements (same y-coordinates)
    V = unique(cols);
    repeated_y_coords = V(histc(cols,unique(cols))>1); % the coordinates
    for cc = 1:length(repeated_y_coords)
        y_inds = find(cluster_elements(:,2)==repeated_y_coords(cc));
        same_col_elements = cluster_elements(y_inds,:);
        same_col_elements = sort(same_col_elements); % sort them according to their x-coordinate
        if size(same_col_elements,1)>1
            for dd = 2:size(same_col_elements)
                dist        = abs(same_col_elements(dd,:)-same_col_elements(dd-1,:));
                distance    = sum(abs(dist));
                if distance == 2
                    curr_y_weight           = (same_col_elements(dd,:)+same_col_elements(dd-1,:))/2;
                    edge_weight_positions   = [edge_weight_positions ; curr_y_weight];
                end
            end
        end
    end
    
    edge_weight_positions = unique(edge_weight_positions,'rows');
    edge_weights = zeros(size(edge_weight_positions,1),1);
    
    for ee = 1:length(edge_weights)
        edge_weights(ee) = weights_and_labels_map(edge_weight_positions(ee,1),edge_weight_positions(ee,2));
    end
    
end

cluster_element_positions = cluster_elements;



end

