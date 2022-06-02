function [est_source_features,glimpse_labels, feature_rates, glimpse_features, glimpse_pair_features, glimpse_pair_sim] = join_similar_glimpses(mean_glimpse_features, feature_mats, params)
%JOIN_SIMILAR_GLIMPSES join glimpses that are similar
% 1) determine how similar the glimpses are to each other, either by
% relative features, relative features+dnn, absolute features+histogram, or
% absolute stacked features+dnn
% 2) source detection by finding larger clusters, either in the histogram
% or in a similarity space
% 3) assign glimpses to sources
% 
% IN 
%       mean_glimpse_features - cell array with averaged features in each
%                               glimpse. Each cell contains a single
%                               feature type (DIMs: number_of_glimpses X
%                               feature_dimensionality e.g. number of
%                               tested azimuths, or f0s, or 1 if this is
%                               the averaged power)
%       feature_mats          - struct with features per single T-F unit
%                               each field contains one a single feature
%                               type. (DIMs: feature_dimensionality X
%                               number_of_time_frames X number_of_subbands)
%       params                - parameter struct
% 
% OUT
%       est_source_features   - features of estimated sources
%       glimpse_labels        - source label for each glimpse (DIMs: number of glimpses X 1)
%       feature_rates         - only if histogram method is used. Rates with which each possible feature value occurs (DIMs: numel(feature_values) X 1)
%       glimpse_features      - average features in each glimpse (DIMs: number of glimpses X number of features)
%       glimpse_pair_features - only if using stacked or relative glimpse features 
%                               holds stacked/relative features for each pair of 2 glimpses(DIMs: number of glimpse pairs X length of feature vector)
%       glimpse_pair_sim      - only if using stacked or relative glimpse features in combination with a dnn
%                               hold a similaritiy value between 0 and 1 for each glimpse pair (DIMs: number of glimpse pairs X 1)


group_by                    = params.global_params.group_by; % 'absolute_glimpse_features','absolute_glimpse_similarity'
grouping_method             = params.global_params.grouping_method; % 'histogram1D','histogram2D','none','linkage','k_means'
feature_types               = params.global_params.feature_types;

%% determine glimpse similarity
if strcmp(group_by,'absolute_glimpse_features')==1   % only absolute features per glimpse, can be used e.g. in combination with histogram grouping    
    glimpse_features        = mean_glimpse_features{1,1};

    glimpse_pair_features   = [];
    glimpse_pair_sim        = [];
    
elseif strcmp(group_by,'stacked_glimpse_features') == 1 % absolute features per glimpse are stacked (for glimpse pairs) and dnn is used to predict similarity of the two glimpses     
    num_feature_types       = size(params.global_params.feature_types,1);
    glimpse_features_norm   = cell(size(mean_glimpse_features));
    for jj = 1:num_feature_types
        normFeatureTypes            = feature_types(jj,:);
        glimpse_features_norm{jj}   = normSingleFeatureDims(mean_glimpse_features{jj},normFeatureTypes,params);   % normalize to positive values and/or range 0-1, crucial for training and inference  
    end
    [stacked_features]      = stack_glimpse_features(glimpse_features_norm,glimpse_features_norm, 'same');
    [glimpse_pair_sim]      = use_glimpse_pair_dnn(stacked_features,params); % converts a vector of absolute stacked features into one similarity values for each glimpse pairs
    
    glimpse_features        = glimpse_features_norm;
    glimpse_pair_features   = stacked_features;
    
elseif strcmp(group_by,'relative_glimpse_features')==1 % make glimpse pairs and calculate relative features that reflect the similarity of both glimpses, e.g. spatial distance between the glimpses
    relative_features       = relative_glimpse_features(mean_glimpse_features,mean_glimpse_features,params,'same');
    
    if strcmp(params.global_params.model_confs.model_type,'dnn') == 1 % combine different relative features with a dnn
        relative_features   = norm_glimpse_grouping_features(relative_features,mean_glimpse_features,params);
        glimpse_pair_sim    = use_glimpse_pair_dnn(relative_features,params);
    elseif strcmp(params.global_params.model_confs.model_type,'none') == 1 % use a single relative feature without dnn, or use this case for generating training data
        relative_features   = norm_glimpse_grouping_features(relative_features,mean_glimpse_features,params);
        glimpse_pair_sim    = relative_features;
    end
    
    glimpse_features        = mean_glimpse_features;
    glimpse_pair_features   = relative_features;
    
else
    error('define feature conversion')
end

%% source detection and glimpse assignment to sources

if strcmp(grouping_method,'histogram1D') == 1  % histogram method, this is always combined with: group_by='absolute_glimpse_features'
    % histogram method requires the actual values which the regarded feature can take = hist_options.feature_vals
    if strcmp(params.global_params.feature_types{1,1},'mean_nac') == 1
        period_vec                      = params.feature_params.minLags : 1: params.feature_params.maxLags;
        f0_vec                          = (1./period_vec)*params.prep_params.fs; % fundamental frequencies in Hz 
        hist_options.feature_vals       = f0_vec;
        hist_options.tolerance_radius   = params.global_params.pitch_tolerance_radius;
        hist_options.hist_peaks         = params.global_params.hist_peaks;      % 'per_t_f_unit','per_time_frame','per_glimpse'
        feats                           = feature_mats.complete_NAC_mat;
    elseif strcmp(params.global_params.feature_types{1,1},'mean_log_azm') == 1 || strcmp(params.global_params.feature_types{1,1},'peak_azm_val') == 1 || strcmp(params.global_params.feature_types{1,1},'reliable_peak_azm') == 1
        hist_options.feature_vals       = params.global_params.azm_vec;         %  possible spatial positions in ° azimuth
        hist_options.tolerance_radius   = params.global_params.azm_tolerance_radius;
        feats                           = feature_mats.log_az_probs;
        hist_options.hist_peaks         = params.global_params.hist_peaks;      % 'per_t_f_unit','per_time_frame','per_glimpse'
    end
    hist_options.est_num_sources        = params.global_params.est_num_sources;    
    % create a histogram to find which feature values occur most frequently, then determine which glimpses are close to these feature values
    [est_source_features, feature_rates, glimpse_labels] = assign_to_hist_peaks1D(feats,glimpse_features,hist_options);

elseif strcmp(grouping_method,'linkage') == 1               % uses hierarchical clustering   
    square_glimpse      = squareform(glimpse_pair_sim);         % similarity of glimpse pairs reshaped for input to linkage
    tree                = linkage(square_glimpse,'complete');   
    glimpse_labels      = cluster(tree,'maxclust',params.global_params.est_num_sources); 
    
    [N,edges]           = histcounts(glimpse_pair_sim,[0:0.01:1]);
    feature_rates       = cat(2,edges,N);
    est_source_features = [];
    
elseif strcmp(grouping_method,'k_means') == 1 % use k-means to find cluster of similar glimpses
    square_glimpse      = squareform(glimpse_pair_sim);
    glimpse_labels      = kmeansSTEV(square_glimpse,params.global_params.est_num_sources);
    
    [N,edges]           = histcounts(glimpse_pair_sim,[0:0.01:1]);
    feature_rates       = cat(2,edges,N);
    est_source_features = [];    
    
elseif strcmp(grouping_method,'k_means_weight') == 1 % use k-means variant to find cluster of similar glimpses
    square_glimpse      = squareform(glimpse_pair_sim);
    glimpse_labels      = kmeansw(square_glimpse,params.global_params.est_num_sources,big_glimpse_sizes);
    
    [N,edges]           = histcounts(glimpse_pair_sim,[0:0.01:1]);
    feature_rates       = cat(2,edges,N);
    est_source_features = [];

elseif strcmp(grouping_method,'k_meansSARI') == 1 % use k-means variant to find cluster of similar glimpses
    square_glimpse      = squareform(glimpse_pair_sim);
    glimpse_labels      = kmeansSARI(square_glimpse,params.global_params.est_num_sources);  %  glimpse_labels = kmeansSTEV(square_glimpse,params.global_params.est_num_sources);
    
    [N,edges]           = histcounts(glimpse_pair_sim,[0:0.01:1]);
    feature_rates       = cat(2,edges,N);
    est_source_features = [];       

elseif strcmp(grouping_method,'none') == 1 % skip the grouping of glimpses, e.g. when generating training data or just doing ROC analysis  
    est_source_features = [];
    glimpse_labels      = []; 
    feature_rates       = []; 

else
    error('define glimpse grouping method')
end


end

