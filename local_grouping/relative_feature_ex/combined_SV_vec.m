function [ SV_mat_comb, SV_horz_comb, SV_vert_comb ] = combined_SV_vec( feature_mats, params, mask_opts, mode )
%COMBINED_SV_VEC Calculate different types of relative features, here
%referred to as SVs (Similarity Values). The SV reflects how similar the features in
%neighboring T-F units units are.
% 
% INPUTS : features_all - structure containing all features extracted per
%                       t-f-unit 
%           params      - structure with all different types of parameter
%                       settings
%           opts        - options (use binary mask y/n,type of binary mask)
%           mode        - defines the types of SV where each SV is the
%                       combination of a FEATURE type and a CORRELATION
%                       method (this can also be the absolute difference or
%                       the sum etc.)
%                       what the input cell (mode) should look like (for example):
%                       mode      = {   'periodicity_degree','pearsons_corr';...
%                                       'normalized_AC'     ,'general_corr';...
%                                       'power'             ,'abs_diff'};
% OUTPUTS: SV_mat_comb  - matrix with SVs for all transitions and
%                       placeholders for TF units and vertices 
%                       dimensions: (num_win*2-1) X (nFilts*2-1) X num_SV_types
%         SV_horz_comb  - all SVs for within frequency band transitions
%         SV_vert_comb  - all SVs for across frequency band transitions
%
%           
% NOTE: The order of SV-types must always remain the same for consistency
% between inference and trained features vectors. To ensure the correct order of
% the feature types, the function 'SVs_control_order' can be used anytime
% 


%% put the cell entries into the right order to guarantee they are matching the pre-trained GMMs

num_SV_types = size(mode,1);


all_fields = fieldnames(feature_mats);
first_field = all_fields{1,1};
long_name = ['feature_mats.',first_field];

if strcmp(first_field,'power_mat')==1  % in this case the dimensions of the feature_mat are (num_win X nFilts)
    num_win = size(eval(long_name),1);
    nFilts = size(eval(long_name),2);
else                                        % in all other cases the dims are (length(feature_vec) X num_win X nFilts)
    num_win = size(eval(long_name),2);
    nFilts = size(eval(long_name),3);
end



SV_mat_comb          = ones((num_win*2-1), (nFilts*2-1),num_SV_types);   % matrix with similarity values (SV) of next neighbours
SV_vert_comb         = zeros(num_win, nFilts-1, num_SV_types);             % matrix with all across frequency band Similarity Values
SV_horz_comb         = zeros(num_win-1, nFilts, num_SV_types);           % matrix with all within frequency band Similarity Values

for aa = 1:num_SV_types
    %************************************************ feature type
    if strcmp(mode{aa,1},'periodicity_degree')==1
        temp_features = feature_mats.complete_PD_mat;
    elseif strcmp(mode{aa,1},'normalized_ac')==1
        temp_features = feature_mats.complete_NAC_mat;
    elseif strcmp(mode{aa,1},'comb_filter_ratio')==1
        temp_features = feature_mats.complete_CFR_mat;
    elseif strcmp(mode{aa,1},'log_azimuth_probabilities')==1
        temp_features = feature_mats.log_az_probs;
    elseif strcmp(mode{aa,1},'azimuth_probabilities')==1
        temp_features = feature_mats.az_probs;
    elseif strcmp(mode{aa,1},'time_sig')==1
        temp_features = feature_mats.time_sig;
    elseif strcmp(mode{aa,1},'power')==1
        temp_features = mean(feature_mats.power_mat,3); % averaged across both ear channels
        % normalize to a range between 0 and 1
        %     features = (features-min(min(features)))./(max(max(features))-min(min(features)));        
    elseif strcmp(mode{aa,1},'log_power')==1
        temp_features = 10*log10(mean(feature_mats.power_mat,3)); % averaged across both ear channels   
    elseif strcmp(mode{aa,1},'croot_power')==1 % cube root compressiong
        temp_features = (mean(feature_mats.power_mat,3)).^(1/3); % averaged across both ear channels   
    else
        error('undefined feature type for local clustering')
    end
    
    %********************************************** check bin mask options
    if strcmp(mask_opts{1},'lc_WITH_bin_mask')==1                                  % 'y' yes apply the provided binary mask
        bin_mask = mask_opts{2};
        temp_features = apply_bin_mask(temp_features, bin_mask);
    elseif strcmp(mask_opts{1},'lc_WITHOUT_bin_mask')==1                              % 'n' no binary mask applied
        bin_mask = ones(num_win,nFilts);
    else
        disp('should a binary mask be applied before local clustering?')
    end
    
    %******************************************* calculate the correlation
    if strcmp(mode{aa,2},'pearsons_corr') == 1
        for dd = 1:(nFilts-1)                               % loop over filterbands (1st calculate all the similarity values across filterbands),same time frame, different filterband
            for jj = 1:num_win                             
                corr1                           = pearsons_corr(temp_features(:,jj,dd),temp_features(:,jj,(dd+1)));
                corr                            = (corr1+1)/2;
                SV_mat_comb(jj*2-1,dd*2,aa)     = corr;
                SV_vert_comb(jj,dd,aa)          = corr;
            end
        end
        for kk = 1:nFilts                                   % loop over all filterbands (2nd calculate all the similarity values within filterbands), same filterband, different time frame
            for ll = 1:(num_win-1)                          
                corr1                           = pearsons_corr(temp_features(:,ll,kk),temp_features(:,(ll+1),kk));  
                corr                            = (corr1+1)/2;
                SV_mat_comb(ll*2,kk*2-1,aa)     = corr;
                SV_horz_comb(ll,kk,aa)          = corr;
            end
        end
        
    elseif strcmp(mode{aa,2},'general_corr') == 1              % like above (pearsons corr), coulod maybe try using xcorr or something (xcorr(temp_features(:,ll,kk),temp_features(:,(ll+1),kk)))
        for dd = 1:(nFilts-1)
            for jj = 1:num_win
                corr1                           = general_corr(temp_features(:,jj,dd),temp_features(:,jj,(dd+1))); % theoretically, these can be values between -1 and 1
                corr                            = (corr1+1)/2; % restrict to the range between 0 and 1
                SV_mat_comb(jj*2-1,dd*2,aa)     = corr;
                SV_vert_comb(jj,dd,aa)          = corr;
            end
        end
        for kk = 1:nFilts
            for ll = 1:(num_win-1)
                corr1                           = general_corr(temp_features(:,ll,kk),temp_features(:,(ll+1),kk));
                corr                            = (corr1+1)/2;
                SV_mat_comb(ll*2,kk*2-1,aa)     = corr;
                SV_horz_comb(ll,kk,aa)          = corr;
            end
        end     
        
    elseif strcmp(mode{aa,2},'mean_max') == 1               % mean value of the max in both units. this could be used for pitch salience for example
        temp_features                           = (temp_features-min(temp_features(:)))./(max(temp_features(:))-min(temp_features(:))); % restrict to a range from 0-1 so that SV is in a range from 0-1 as well        
        for dd = 1:(nFilts-1)                               %  across filterbands
            for jj = 1:num_win                              
                mean_max                        = mean([max(temp_features(:,jj,dd)),max(temp_features(:,jj,(dd+1)))]);              
                SV_mat_comb(jj*2-1,dd*2,aa)     = mean_max;
                SV_vert_comb(jj,dd,aa)          = mean_max;
            end
        end
        for kk = 1:nFilts                                   % within filterbands
            for ll = 1:(num_win-1)                         
                mean_max                        = mean([max(temp_features(:,ll,kk)),max(temp_features(:,(ll+1),kk))]);
                SV_mat_comb(ll*2,kk*2-1,aa)     = mean_max;
                SV_horz_comb(ll,kk,aa)          = mean_max;
            end
        end       
        
    elseif strcmp(mode{aa,2},'abs_diff') == 1               % used for example for the power deltas       
        abs_diffs_vert                          = zeros(num_win,nFilts-1);% for this SV the dimensions must be (num_time_frames X num_filterbands)
        for dd = 1:(nFilts-1)                                   % across filterbands
            for jj = 1:num_win
                abs_diffs_vert(jj,dd)           = abs(temp_features(jj,dd+1)-temp_features(jj,dd)); 
            end
        end
        abs_diffs_horz                          = zeros((num_win-1),nFilts);
        for kk = 1:nFilts                                   % within band
            for ll = 1:(num_win-1)
                abs_diffs_horz(ll,kk)           = abs(temp_features(ll+1,kk)-temp_features(ll,kk));
            end
        end
        
        % norm options
        allover_max = max(max(abs_diffs_vert(:)),max(abs_diffs_horz(:)));
        abs_diffs_vert = abs_diffs_vert./allover_max;
        abs_diffs_horz = abs_diffs_horz./allover_max;

        for dd = 1:(nFilts-1)                 
            for jj = 1:num_win
                SV_vert_comb(jj,dd,aa)          = 1-abs_diffs_vert(jj,dd); % by subtracting the difference from 1, the maximal similarity (=0 difference) yields and SV of 1
                SV_mat_comb(jj*2-1,dd*2,aa)     = 1-abs_diffs_vert(jj,dd);
            end
        end        

        for kk = 1:nFilts 
            for ll = 1:(num_win-1)
                SV_horz_comb(ll,kk,aa)          = 1-abs_diffs_horz(ll,kk);
                SV_mat_comb(ll*2, kk*2-1, aa)   = 1-abs_diffs_horz(ll,kk);
            end
        end
        
    elseif strcmp(mode{aa,2},'abs_sum') == 1
        abs_sums_vert                           = zeros(num_win,nFilts-1); % for this SV the dimensions must be (num_time_frames X num_filterbands)
        for dd = 1:(nFilts-1)                               % across filterbands
            for jj = 1:num_win
                abs_sums_vert(jj,dd)            = abs(temp_features(jj,dd+1)+temp_features(jj,dd));

            end
        end

        abs_sums_horz                           = zeros((num_win-1),nFilts);
        for kk = 1:nFilts                                               % within band
            for ll = 1:(num_win-1)
                abs_sums_horz(ll,kk)            = abs(temp_features(ll+1,kk)+temp_features(ll,kk));

            end
        end
        
        
        % norm options
        allover_max = max(max(abs_sums_vert(:)),max(abs_sums_horz(:)));
        abs_sums_vert = abs_sums_vert./allover_max;
        abs_sums_horz = abs_sums_horz./allover_max;


        for dd = 1:(nFilts-1)                                       % across filterbands
            for jj = 1:num_win
                SV_vert_comb(jj,dd,aa)          = 1-abs_sums_vert(jj,dd); % by subtracting the difference from 1, the maximal similarity (=0 difference) yields and SV of 1
                SV_mat_comb(jj*2-1,dd*2,aa)     = 1-abs_sums_vert(jj,dd);
            end
        end    
        for kk = 1:nFilts  % within band
            for ll = 1:(num_win-1)
                SV_horz_comb(ll,kk,aa)          = 1-abs_sums_horz(ll,kk);
                SV_mat_comb(ll*2, kk*2-1, aa)   = 1-abs_sums_horz(ll,kk);
            end
        end
        
    else
        disp('unknown correlation mode for local clustering')
        
    end
    
        
    
    
end






end

