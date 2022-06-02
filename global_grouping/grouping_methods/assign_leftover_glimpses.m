function [small_glimpse_labels] = assign_leftover_glimpses(small_glimpse_features,big_glimpse_features,big_glimpse_labels,est_source_features,params)
%ASSIGN_LEFTOVER_GLIMPSES if any glimpses are leftover (e.g. because they
%are very small or contain unreliable features etc.), they are assigned to
%the already detected sources. This can either be done by assigning them
%the same source as a larger (already assigned) glimpse, to which they are
%very similar, or by comparing them to the averaged features of the
%estimated sources (est_source_features) and assigning the small glimpse to
%the most similar source.
%   INPUTS:
%       small_glimpse_features  - averaged features for each of the
%                                 leftover (small) glimpses
%       big_glimpse_features    - averaged features for the already
%                                 assigned (big) glimpses
%       est_source_features     - average features of the estimated sources
%       big_glimpse_labels      - labels of the already assigned glimpses
%       params                  - parameter struct
%
%   OUTPUTS:
%       small_glimpse_labels    - new labels for the leftover (small) glimpses


num_big_glimpses    = size(big_glimpse_features{1,1},1);
num_small_glimpses  = size(small_glimpse_features{1,1},1);

%% three different methods to label leftover glimpses, either

if strcmp(params.global_params.group_by,'absolute_glimpse_features')==1 % histogram method
%     [small_abs_glimpse_features1]    = absolute_glimpse_features(small_glimpse_features,params.global_params.feature_types,params);
    small_abs_glimpse_features      = small_glimpse_features{1};
    
    small_glimpse_labels            = zeros(num_small_glimpses,1);
    
    if strcmp(params.global_params.feature_types{1,1},'mean_nac') == 1
        tolerance_radius   = params.global_params.pitch_tolerance_radius;   % tolerance radius around the histogram peak, in which the glimpse is group to the peak
    elseif strcmp(params.global_params.feature_types{1,1},'mean_log_azm') == 1 || strcmp(params.global_params.feature_types{1,1},'peak_azm_val') == 1
        tolerance_radius   = params.global_params.azm_tolerance_radius;
    elseif strcmp(params.global_params.feature_types{1,1},'spectral_centroid') == 1
        tolerance_radius   = 1;
    end
    
    for pp = 1:num_small_glimpses                                              % for each glimpse: check if it is near an overall_feature_peak
        for qq = 1:(params.global_params.est_num_sources-1)                                  % if so: assign the glimpse to the peaks' location
            curr_peak     = est_source_features(qq);
            %curr_peak_angle     = overall_feature_peaks(end-(qq-1));   % 4.2. try this to start with the lower peaks, this gives the higher peaks an advantage
            if (curr_peak-tolerance_radius) <= small_abs_glimpse_features(pp) && small_abs_glimpse_features(pp) <= curr_peak+tolerance_radius
                small_glimpse_labels(pp)      = qq;                           % assign the glimpse an ID that now corresponds to this location
            else
                continue                                                % glimpse_ID_vec(pp) remains zero and is assigned to the 'leftovers' class
            end
        end
        if small_glimpse_labels(pp)==0                                        % if this glimpse was not yet assigned any location
            small_glimpse_labels(pp)      = params.global_params.est_num_sources;                  % is assigned to the 'trash class'(leftovers)
        end
    end

elseif strcmp(params.global_params.group_by,'stacked_glimpse_features')==1   % determine similarity based on stacked feature vectors + dnn
    feature_types               = params.global_params.feature_types;
    num_feature_types       = size(params.global_params.feature_types,1);
    small_glimpse_features_norm   = cell(size(small_glimpse_features));
    big_glimpse_features_norm   = cell(size(big_glimpse_features));
    for jj = 1:num_feature_types
        normFeatureTypes                = feature_types(jj,:);
        small_glimpse_features_norm{jj}   = normSingleFeatureDims(small_glimpse_features{jj},normFeatureTypes,params);
        big_glimpse_features_norm{jj}   = normSingleFeatureDims(big_glimpse_features{jj},normFeatureTypes,params);
    end
    [stacked_features] = stack_glimpse_features(small_glimpse_features_norm,big_glimpse_features_norm, 'different');
    [smallbig_pair_SVs] = use_glimpse_pair_dnn(stacked_features,params);
    
    
    %%%%%%
    
    
%     [small_abs_glimpse_features]        = absolute_glimpse_features(small_glimpse_features,params.global_params.feature_types,params);
%     [big_abs_glimpse_features]          = absolute_glimpse_features(big_glimpse_features,params.global_params.feature_types,params);
%   
%     small_abs_glimpse_features          = norm_glimpse_grouping_features(small_abs_glimpse_features,small_glimpse_features,params); % normalize to positive values and/or range 0-1
%     big_abs_glimpse_features            = norm_glimpse_grouping_features(big_abs_glimpse_features,small_glimpse_features,params); % normalize to positive values and/or range 0-1
%     
%     [smallbig_pair_SVs, ~]              = stacked_features_sim(small_abs_glimpse_features,big_abs_glimpse_features,params,'different'); % smallbig_pair_SVs has the dimensions (small glimpsesXbig glimpses)
%     
    smallbig_mat                        = reshape(smallbig_pair_SVs,num_big_glimpses,num_small_glimpses); % similarities between the small and big glimpses
    [~,small_labels_inds]               = max(smallbig_mat,[],1); % determine to which of the big glimpses each small glimpse is most similar (max value in the smallbig_pair_SVs matrix)
    small_glimpse_labels                = big_glimpse_labels(small_labels_inds);
    
    
    
elseif strcmp(params.global_params.group_by,'relative_glimpse_features')==1   % determine similarity based on relative features and optionally a dnn
    [relative_features]         = relative_glimpse_features(small_glimpse_features,big_glimpse_features,params,'different');
    
    if strcmp(params.global_params.model_confs.model_type,'dnn') == 1 % possibility to combine different relative features with a dnn
        relative_features                       = norm_glimpse_grouping_features(relative_features,small_glimpse_features,params);
        glimpse_pair_sim                        = glimpse_sim_relative_features(relative_features,params);
    elseif strcmp(params.global_params.model_confs.model_type,'none') == 1 % use a single relative feature without dnn, or use this case for generating training data
        relative_features                       = norm_glimpse_grouping_features(relative_features,small_glimpse_features,params);
        
        glimpse_pair_sim                        = relative_features;
    end
    
    smallbig_mat                    = reshape(glimpse_pair_sim,num_big_glimpses,num_small_glimpses); % similarities between the small and big glimpses
    [~,small_labels_inds]           = max(smallbig_mat,[],1); % determine to which of the big glimpses each small glimpse is most similar (max value in the smallbig_pair_SVs matrix)
    small_glimpse_labels            = big_glimpse_labels(small_labels_inds);

else
    error('define feature conversion')
end




end

