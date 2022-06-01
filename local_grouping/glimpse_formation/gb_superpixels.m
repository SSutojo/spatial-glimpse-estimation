function [contours, glimpse_map, weights_and_labels_map] = gb_superpixels(contrast_mat, tau_param)
%USE A GRAPH-BASED SUPERPIXELS TO DIVIDE IMAGE INTO REGIONS WITH SIMILAR FEATURES
% Felzenszwalb, P., Huttenlocher, D.P.(2004),"Efficient Graph-Based Image
% Segmentation", International Journal of Computer Vision 59(2),167-181
% 
% This function gives every t-f-unit in weight_mat its own label, then goes
% through all edge weights (sorted_weights2), finds the neighboring t-f-units
% and gets their labels. Then calculates the inner distance for the two clusters
% compares the inner distance with the regarded edge weight and merges or
% doesn't merge by giving the two regions a new label (start from 1, set a counter or something )
% 
% 
% INPUTS:
%   contrast_mat    - matrix that contains all contrast values betwenn neighboring pixels. 
%                     Dimensions: (time_frames*2-1 , frequency_bands*2-1)
%                     zeros in places where a t-f-unit is and edge weights between them
%   tau_param       - parameter to control the size of the formed superpixels
%
% OUTPUTS:
%   contours        - matrix of same size and structure as contrast_map but
%                     with binary values, indicating the contours of
%                     superpixels
%   glimpse_map     - 2dim matrix (number of time frames X number of filters)
%                     each TF unit is labeled with an integer, if two TF
%                     units are labeled with the same integer, they belong
%                     to the same superpixel
%   weights_and_labels_map - same dimensions as contrast_mat, but this
%                     contains also the integers indicating the formed
%                     superpixels
%



weight_mat       = contrast_mat;            % weight_mat = 1-SV_mat; % the larger the weight, the lower the similarity,  the colums are joint (1st col, 2nd col, 3rd col...)

num_bands       = ceil(size(contrast_mat,2)/2);
num_frames      = ceil(size(contrast_mat,1)/2);
num_pixels      = num_bands * num_frames;           % equal to the number of t-f-units


%%  sort edge weights, starting with the smallest weight (=highest similarity ) and remove placeholders for t-f-units and vertexes
% because in weight_mat only every other element (even number) contains an actual weight and every odd element is a placeholder for 
% a t-f-unit (or a vertex with no specified value) every odd element from sorted_weight is removed

all_weights_vec                             = weight_mat(:);   % edge_weights vector,
[sorted_weights, sorted_weights_indices]    = sort(all_weights_vec,'ascend');  

odd_element_indices                         = find(mod(sorted_weights_indices,2)); % for odd numbers  modulus is 1                                                     
sorted_weights2                             = sorted_weights;
sorted_weights_indices2                     = sorted_weights_indices;

% remove odd indices from  sorted_weights and sorted weight_indices
sorted_weights2(odd_element_indices)        = [];  % these are the weights that should be considered
sorted_weights_indices2(odd_element_indices)= [];  % these are the indices corresponding to the considered weights (in weight_mat)


%% initialize new contour map (weights_and_labels_map)
label_vec               = (num_pixels+1):2*num_pixels;              % start with every t-f-unit having an individual label (=temporary labels)
label_mat               = reshape(label_vec,num_frames,num_bands);  % these are initial labels to be inserted into initial weights_and_labels_map

% positions of t-f-units in the weights_and_labels_map
frames_vec              = 1:2:(2*num_frames-1);
bands_vec               = 1:2:(2*num_bands-1);                     

weights_and_labels_map  = weight_mat;                               % weights taken from weight_mat, every t-f-unit gets a new "superpixel label"
weights_and_labels_map(frames_vec,bands_vec) = label_mat;           % insert the inital labels


%% starting with the lowest weight, compare neighboring clusters (at this edge weight) and decide wether to fuse them or draw boundary
for aa=1:length(sorted_weights2)
    % choose on edge weight and find the current position of the regarded edge in weights_and_labels_map 
    % by converting the index into a frequency band and frame position
    current_weight                  = sorted_weights2(aa);
    current_index                   = sorted_weights_indices2(aa);    
    [current_row, current_column]   = ind2sub(size(weights_and_labels_map),current_index);
    
    %% compare the neighboring clusters at the current edge weight
    % if current column is an odd number, look on top and below (+/- one time frame, 1st dimension)
    % if current column is an even number, look on left and right (+/- one frequency band, 2nd dimension)    
    if mod(current_column,2)==1                                             % odd number
        % get the labels of the two pixels connected by the regarded edge
        label_A = weights_and_labels_map(current_row-1, current_column);
        label_B = weights_and_labels_map(current_row+1, current_column);        
        if label_A ~= label_B                                               % continue if the pixels have different labels
            % get all weights and numbers of elements in both clusters
            [edge_weights_A, cluster_size_A , cluster_element_positions_A] = get_edge_weights(weights_and_labels_map, label_A);
            [edge_weights_B, cluster_size_B , cluster_element_positions_B] = get_edge_weights(weights_and_labels_map, label_B);
            
%             if cluster_size_A > max_size || cluster_size_B>= max_size
%                 continue
%             end
            % calculate the minimum internal difference "MInt" and compare to current weight
            tau_A       = tau_param/cluster_size_A;
            int_A       = max(edge_weights_A) + tau_A;
            tau_B       = tau_param/cluster_size_B;
            int_B       = max(edge_weights_B) + tau_B;
            MInt        = min(int_A,int_B);
            if current_weight > MInt                    % create a boundary, do not merge the clusters               
                continue
            elseif current_weight <= MInt               % merge the two clusters                               
                % choose the lower integer (as new label) to overwrite the higher integer and write the same labels into both clusters
                % turn cluster_elements_positions into linear indices of weights_and_labels_map
                new_label       = min(label_A, label_B);
                linear_inds_A   = sub2ind(size(weights_and_labels_map),cluster_element_positions_A(:,1), cluster_element_positions_A(:,2));                
                linear_inds_B   = sub2ind(size(weights_and_labels_map),cluster_element_positions_B(:,1), cluster_element_positions_B(:,2));
                weights_and_labels_map(linear_inds_A) = new_label;
                weights_and_labels_map(linear_inds_B) = new_label;               
            end
        end
        
        
    elseif mod(current_column,2)~=1                                         % even number
        label_A = weights_and_labels_map(current_row, current_column-1);
        label_B = weights_and_labels_map(current_row, current_column+1);
        if label_A~=label_B
            [edge_weights_A, cluster_size_A , cluster_element_positions_A] = get_edge_weights(weights_and_labels_map, label_A );
            [edge_weights_B, cluster_size_B , cluster_element_positions_B] = get_edge_weights(weights_and_labels_map, label_B);
%             if cluster_size_A > max_size || cluster_size_B>= max_size
%                 continue
%             end
            
            tau_A   = tau_param/cluster_size_A;
            int_A   = max(edge_weights_A) + tau_A;
            tau_B   = tau_param/cluster_size_B;
            int_B   = max(edge_weights_B) + tau_B;
            MInt    = min(int_A,int_B);
            if current_weight > MInt
                continue
            elseif current_weight <= MInt
                new_label       = min(label_A, label_B);
                linear_inds_A   = sub2ind(size(weights_and_labels_map),cluster_element_positions_A(:,1), cluster_element_positions_A(:,2));
                linear_inds_B   = sub2ind(size(weights_and_labels_map),cluster_element_positions_B(:,1), cluster_element_positions_B(:,2));
                weights_and_labels_map(linear_inds_A)= new_label;
                weights_and_labels_map(linear_inds_B)=new_label; 
            end
        end
         
    end   
end 

%% convert weights_and_labels_map into glimpse_map and contour map
glimpse_map         = zeros(num_frames, num_bands);
for bb = 1:num_frames
    for cc = 1:num_bands
        glimpse_map(bb,cc) = weights_and_labels_map((2*bb)-1,(2*cc)-1);
    end
end

% glimpses wieder von 1 beginnend durchnummerieren
old_glimpse_IDs     = unique(glimpse_map);
old_glimpse_IDs     = sort(old_glimpse_IDs);
new_glimpse_matrix  = zeros(size(glimpse_map,1), size(glimpse_map,2),length(old_glimpse_IDs));
for cc = 1:length(old_glimpse_IDs)
    new_glimpse_matrix(:,:,cc) = (glimpse_map==old_glimpse_IDs(cc)).*cc;
end
new_glimpse_map     = sum(new_glimpse_matrix,3);
glimpse_map         = new_glimpse_map;

% get the contours of the glimpse map
contours = get_contours(glimpse_map);


end