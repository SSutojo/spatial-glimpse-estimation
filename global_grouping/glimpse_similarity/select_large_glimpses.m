function [num_big_glimpses, num_small_glimpses, big_glimpse_features, small_glimpse_features, big_swapIDs, small_swapIDs , big_glimpse_sizes] = select_large_glimpses(glimpse_sizes,min_glimpse_size, mean_glimpse_features)
%SELECTLARGEGLIMPSES from a vector of glimpses with different sizes,
%selects the large enough glimpses ('size' refers to the number of T-F
%units in the glimpse)
%
% INPUTS 
%   glimpse_size            : vector (length=number of glimpses) which contains the size of each glimpse
%   min_glimpse_size        : the minimum glimpse size (smaller glimpses are excluded)
%   mean_glimpse_features   : averaged features for each glimpse, cell array
%                             dims: number of feature types X 1
%                             each cell contains 1 feature type (num. of glimpses X length of feature vector)
% 
% OUTPUTS
%   num_big_glimpses        : number of glimpses above the minimum glimpse size
%   num_small_glimpses      : number of leftover (smaller) glimpses
%   big_glimpse_features    : cell array with features of big glimpses (similar dimensions as mean_glimpse_features)
%   small_glimpse_features  : cell array with features of small glimpses (similar dimensions as mean_glimpse_features)
%   big_swapIDs             : 2 column vector (left column: old glimpse ID (relevant to original set of all glimpses); right column: new glimpse ID)
%   small_swap_IDs          : 2 column vector (left column: old glimpse ID (relevant to original set of all glimpses); right column: new glimpse ID)
%   big_glimpse_sizes       : sizes of big glimpses

%% find all glimpses above minimum size
big_glimpse_IDs             = find(glimpse_sizes>=min_glimpse_size);  % z.B.mindest glimpse size ist 5 T-F units
big_swapIDs                 = [big_glimpse_IDs,[1:length(big_glimpse_IDs)].'];  % 1st column contains original IDs, 2nd column the new IDs
small_glimpse_IDs           = find(glimpse_sizes < min_glimpse_size);
small_swapIDs               = [small_glimpse_IDs,[1:length(small_glimpse_IDs)].'];

%% for each feature type, select only the big glimpses and collect in a new cell array
num_feature_types           = size(mean_glimpse_features,1);
big_glimpse_features        = cell(num_feature_types,1);
small_glimpse_features      = cell(num_feature_types,1);

for aa = 1:num_feature_types  
    curr_glimpse_feats          = mean_glimpse_features{aa};
    big_glimpse_features{aa}    = curr_glimpse_feats(big_glimpse_IDs,:);
    small_glimpse_features{aa}  = curr_glimpse_feats(small_glimpse_IDs,:);
    
end
num_big_glimpses            = numel(big_glimpse_IDs);
num_small_glimpses          = numel(small_glimpse_IDs);
big_glimpse_sizes           = glimpse_sizes(big_glimpse_IDs);


end

