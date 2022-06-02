function [main_out] = spatial_glimpses_main(noisy_mix,params)
%SPATIAL_GLIMPSES_MAIN extracts spatial glimpses and estimates binary masks
%for signal noisy_mix. For more details see subfunctions.
% 
% IN
%   noisy_mix   - binaural input signal 
%   params      - parameter struct 
% OUT
%   main_out    - results struct with features, estimates masks, etc.


% ******************* pre-processing (division into time frames and filterbands) *********************************
[pre_processed]             = pre_process(noisy_mix, params.prep_params);

% ******************* extract features  **************************************************************************
[feature_mats]              = feature_ex(pre_processed, {'PERIODICITY','AZIMUTH','SUMMED_POWER'}, params.feature_params, params.prep_params);

% ******************** local clustering (glimpse formation) ******************************************************
[clustering_out]            = local_clustering(feature_mats,params);

% ******************** global clustering (glimpse labeling) *************
[global_out]                = global_grouping(clustering_out.glimpse_map, feature_mats, params);


main_out.pre_processed  = pre_processed;
main_out.feature_mats   = feature_mats;
main_out.clustering_out = clustering_out;
main_out.global_out     = global_out;
end

