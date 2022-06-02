function [ out_struct ] = global_grouping( glimpse_map, feature_mats, params )
%GLOBAL_GROUPING assign the local clusters to sources and get an estimate of the source map 
%
% - average features within each glimpse and determine glimpse size
% - join similar glimpses, either by using glimpse similarities (can be based on relative fetures, or stacked absolute features + a DNN) 
%   or by creating a histogramm with feature rates and assigning glimpses to the most frequently appearing feature values
% - create an estimated source map by inserting estimated glimpse IDs into
%   estimated glimpse map
%   
% IN:   
%   glimpse_map - map in which every time-frequency unit is labeled with an integer that indicates
%                 which glimpse it belongs to (DIMs: num. of time frames X num. of subbands) 
%   feature_mats- structure that holds different features types/ feature matrices.
%                 Each feature matrix has dimensions (length_of_feature_vector X num_win X num_filterbands)
%   params      - parameter struct
%
% OUT: 
%   out_struct.  ... see more detailed description at the end of this function
%               

num_glimpses                        = max(glimpse_map(:));

%% average features over glimpses
[feature_types]                     = glob_control_order(params.global_params.feature_types);% make sure the features are always in the same order
params.global_params.feature_types  = feature_types;

% optional pre-processing depending on feature type, e.g. emphasize the lower f0 values
if params.global_params.low_f0_emphasis==true
    [feature_mats] = pre_process_glimpse_features(feature_mats,params);   
end

[mean_glimpse_features]             = average_glimpse_features(feature_mats,glimpse_map,feature_types(:,1),params);

%% option to select large glimpses only and for glimpse pairs and clustering
[glimpse_size_cell]                 = average_glimpse_features(feature_mats,glimpse_map,{'glimpse_size'},params);
glimpse_sizes                       = glimpse_size_cell{1,1};
if params.global_params.min_glimpse_size >= 2   % just use glimpses larger than the minimum size
    [~, ~, big_glimpse_features, small_glimpse_features, big_swapIDs, small_swapIDs , ~] = select_large_glimpses(glimpse_sizes,params.global_params.min_glimpse_size, mean_glimpse_features);   
else
    big_glimpse_features = mean_glimpse_features;  % in this case all glimpses are used for the grouping
end

%% determine the similarity of glimpses to each other and group them 
[est_source_features,glimpse_labels, feature_rates, glimpse_features, glimpse_pair_features, glimpse_pair_sim] = join_similar_glimpses(big_glimpse_features, feature_mats, params);
big_glimpse_labels              = glimpse_labels;

%% also assign small glimpses to the detected sources
if params.global_params.min_glimpse_size >= 2
    [small_glimpse_labels]      = assign_leftover_glimpses(small_glimpse_features,big_glimpse_features,big_glimpse_labels,est_source_features,params); 
    % combine small_glimpse_labels and big glimpse labels
    newT                        = zeros(num_glimpses,1);
    newT(big_swapIDs(:,1),1)    = glimpse_labels;   % go back to the before IDs
    newT(small_swapIDs(:,1),1)  = small_glimpse_labels;   
    glimpse_labels              = newT;  
end

%% write glimpse labels into est_source_map 
if strcmp(params.global_params.grouping_method,'none') == 1
    est_source_map = [];
else
    est_source_map                  = zeros(size(glimpse_map));
    for zz = 1:num_glimpses
        [rowz,colz]                 = find(glimpse_map == zz);              % get row/column indices of the elements in the regarded glimpse
        size_zz                     = length(rowz);                         % number of t-f-units units in the regarded glimpse
        for uu = 1:size_zz
            est_source_map(rowz(uu),colz(uu)) = glimpse_labels(zz);
        end
    end
end

%%
out_struct.est_source_map           = est_source_map;           % map in which each T-F unit is labeled with the (estimated) dominant source
out_struct.all_glimpse_features     = mean_glimpse_features;    % average features in each glimpse (DIMs: number of glimpses X number of features)
out_struct.big_glimpse_features     = big_glimpse_features;     % average features of glimpses that are above the minimum glimpse size
out_struct.glimpse_features         = glimpse_features;         % either the same as mean_glimpse_features, or (if group_by='stacked_glimpse_features') normalized and stacked features

% only if using stacked or relative glimpse features:
out_struct.glimpse_pair_sim         = glimpse_pair_sim;         % similaritiy value between 0 and 1 for each glimpse pair (DIMs: number of glimpse pairs X 1)
out_struct.glimpse_pair_features    = glimpse_pair_features ;   % stacked/relative features for each pair of 2 glimpses(DIMs: number of glimpse pairs X length of feature vector)

% only if using histogram method
out_struct.peak_features            = est_source_features;      % features of estimated sources
out_struct.feature_rates            = feature_rates;            % Rates with which each possible feature value occurs (DIMs: numel(feature_values) X 1)


end







