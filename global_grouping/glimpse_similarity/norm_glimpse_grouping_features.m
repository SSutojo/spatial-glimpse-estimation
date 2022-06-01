function [featuresOUT] = norm_glimpse_grouping_features(featuresIN,mean_glimpse_features,params)
%NORM_GLIMPSE__GROUPING_FEATURES Normalize glimpse similarity features for global
%grouping (used in training and inference)
% mean glimpse features: needed to determine feature dimensions

% for relative features
% DIMS INPUT:
% num_datapoints X num_feature_types (each relative feature is dim =1)

% for absolute features
% DIMS INPUT:
% num_datapoints X num_feature_types/feat_lengths (an absolute feature can have dim > 1 )



% determine the size of the 

feature_types = params.global_params.feature_types;

num_feature_types = size(params.global_params.feature_types,1);
feat_lengths                = zeros(num_feature_types,1);
start_idx                   = zeros(num_feature_types,1);
end_idx                     = zeros(num_feature_types,1);

for bb = 1:num_feature_types                                % prepare length of the vector that will later contain all features 
    if strcmp(feature_types{bb,2},'direct') == 1        % the entire feature vector of a glimpse is used (e.g. probabilities of all possible azimuth values)
        feat_lengths(bb)    = size(mean_glimpse_features{bb},2);
    elseif strcmp(feature_types{bb,2},'direct') == 0    % only one feature value is used (e.g. only the azimuth value with maximum probability)
        feat_lengths(bb)    = 1;
    end    
    end_idx(bb)             = sum(feat_lengths(1:bb));
    start_idx(bb)           = end_idx(bb) -(feat_lengths(bb)-1);
end

num_glimpses                = size(mean_glimpse_features{1},1);

featuresOUT       = zeros(size(featuresIN));  % matrix to hold all features (sum(feat_lengths)) for all glimpses   :zeros(num_glimpses,sum(feat_lengths));  

% absolute features (prior to stacking - these are for single glimpses)

for aa = 1:num_feature_types
    
    normFeatureTypes  = feature_types(aa,:);
    
    featuresOUT(:,start_idx(aa):end_idx(aa)) = normSingleFeatureDims(featuresIN(:,start_idx(aa):end_idx(aa)),normFeatureTypes,params);
    
end



% relative features (these are for glimpse pairs)


end

