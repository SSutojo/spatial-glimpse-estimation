function [feature_mats] = pre_process_glimpse_features(feature_mats,params)
%PRE_PROCESS_GLIMPSE_FEATURES Summary of this function goes here
%   Detailed explanation goes here
%
% f0 smoothing
% low f0 emphasis
% 
if params.global_params.low_f0_emphasis==true    % emphasize the lower f0 values
    period_vec                      = params.feature_params.minLags:1:params.feature_params.maxLags;
    params.global_params.f0_vec     = params.prep_params.fs./period_vec;
    
    
    % create a vector to emphasize the lower f0s before averaging them
    [~, min_dist_ind]               = min(abs(params.global_params.f0_vec - params.global_params.lowpass_fc)); % find index that is closest to the cutoff frequency (slope changes from 1 to sqrt)
    emph1                           = 0:(min_dist_ind-1); % first part with linear slope
    change_ind                      = round((length(params.global_params.f0_vec)-min_dist_ind)/2)+min_dist_ind; % index where the slope changes from sqrt function to 0
    emph2                           = ((0:(change_ind-min_dist_ind)).^(1/3)+min_dist_ind); % 2nd part where slope is sqrt function (or 3rd root maybe)
    emph3                           = emph2(end).*ones((length(params.global_params.f0_vec)-(change_ind+1)),1); % 3rd part, slope is 0
    emph_vec                        = [emph1(:);emph2(:);emph3(:)];
    emph_vec                        = emph_vec./max(emph_vec);
    for aa = 1:size(feature_mats.complete_NAC_mat,2)
        for bb = 1:size(feature_mats.complete_NAC_mat,3)
            feature_mats.complete_NAC_mat(:,aa,bb) = feature_mats.complete_NAC_mat(:,aa,bb).*emph_vec;
        end
    end    
end


%% OPTION f0 smoothing  
% %     
% %     options.smooth_f0rates      = params.global_params.smooth_f0rates;
% %     options.smooth_f0s_window   = params.global_params.smooth_f0s_window;
% %     options.pitchsal_weight     = params.global_params.pitchsal_weight;
% %     options.lowpass_fc          = params.global_params.lowpass_fc;  % 400
% %     period_vec                  = params.feature_params.minLags : 1: params.feature_params.maxLags; 
% %     options.f0_vec              = (1./period_vec)*params.prep_params.fs; % fundamental frequencies in Hz
% %     options.est_num_sources     = est_num_sources;
% %     options.tolerance_radius    = params.global_params.pitch_tolerance_radius;
% %     options.low_f0_emphasis     = params.global_params.low_f0_emphasis;
% %     features                    = feature_mats.complete_NAC_mat; % just an example, this could also be PD or CFR
% %     [est_source_map, peak_pitches, glimpse_pitches, f0_rates ] = periodicity_global(glimpse_map, features, options);
% %     




end

