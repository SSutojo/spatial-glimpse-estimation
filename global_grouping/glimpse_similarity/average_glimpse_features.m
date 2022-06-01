function [mean_glimpse_features] = average_glimpse_features(feature_mats,glimpse_map,feature_types,params)
%GET_GLIMPSE_FEATURES average the features in each glimpse (features
%within a glimpse may be assumed to be homogeneous or slowly
%changing over time)
%Note: ensure the right order of feature types before using this function
%(use glob_control_order)
% 
% IN:
%   features_mats   - features in each T-F unit (each cell contains a single feature type)
%   dimensions of each individual feature_mat (dims:length_of_feature_vec X num.frames X num.bands)
%   glimpse_map     - 2 dimensional matrix in which each T-F units
%                     belonging to the same glimpse are labeled with the
%                     same integer
%   feature_types   - the feature types that should be averaged
%   params          - paramter struct
% OUT:
%   mean_glimpse_features - cell array that contains the averaged features
%                     in each glimpse. Dims: number of glimpses X lenght of averaged feature vector
%                     each cell contains a single feature type

num_glimpses                    = max(glimpse_map(:));                      % get number of separate glimpses /clusters of T-F units
num_feature_types               = size(feature_types,1);
mean_glimpse_features           = cell(num_feature_types,1);                % first column is glimpse ID 

for aa = 1:num_feature_types
    if strcmp(feature_types{aa},'mean_log_azm')==1                          % average spatial location, save the vector with all azimuth probabilities
        feature_mat                 = feature_mats.log_az_probs;
        azm_vec                     = params.global_params.azm_vec;
        num_azms                    = numel(azm_vec);
        glimpse_azm_probs           = zeros(num_glimpses,num_azms);
        for zz = 1:num_glimpses
            [rowz,colz]             = find(glimpse_map == zz);              % get row/column indices of the elements in the regarded glimpse
            size_zz                 = length(rowz);                         % number of t-f-units units in the regarded glimpse
            acc_log_azs             = zeros(size_zz,num_azms);              % accumulate all log_like_az vectors beloging to this glimpse in one matrix
            for uu = 1:size_zz
                acc_log_azs(uu,:)   = feature_mat(:,rowz(uu),colz(uu));
            end
            glimpse_azm_probs(zz,:)  = mean(acc_log_azs,1);                 % average over 1st dim  prob. azm 
        end
        mean_glimpse_features{aa} = glimpse_azm_probs;                      % num_glimpses X num_azimuths
        
    elseif strcmp(feature_types{aa},'peak_azm_val')==1                      % average spatial location, save the azimuth value with the highest probability
        feature_mat                 = feature_mats.log_az_probs;
        azm_vec                     = params.global_params.azm_vec;
        num_azms                    = numel(azm_vec);
        peak_azm_vals               = zeros(num_glimpses,1);
        for zz = 1:num_glimpses
            [rowz,colz]             = find(glimpse_map == zz);              % get row/column indices of the elements in the regarded glimpse
            size_zz                 = length(rowz);                         % number of t-f-units units in the regarded glimpse
            acc_log_azs             = zeros(size_zz,num_azms);              % accumulate all log_like_az vectors beloging to this glimpse in one matrix
            for uu = 1:size_zz
                acc_log_azs(uu,:)   = feature_mat(:,rowz(uu),colz(uu));
            end
            mean_probs  = mean(acc_log_azs,1);                              % azimuth probabilities
            [~,Mc,~,~]           = findLocalPeaks(mean_probs.');            % get index of the peak azimuth for each glimpse
            if isempty(Mc)
                [~,Mc]           = max(mean_probs);                         % use the max value if no peak exists
            end
            peak_azm_vals(zz)  =  azm_vec(Mc);    
        end
        mean_glimpse_features{aa} = peak_azm_vals;                          % num_glimpses X 1
        
     elseif strcmp(feature_types{aa},'reliable_peak_azm')==1       % average spatial location, use cross correlation as reliability criterion, save the azimuth value with the highest probability
        feature_mat                 = feature_mats.log_az_probs;
        maxXCmat                    = max(feature_mats.xcorrs,[],1); % max. of cross correlation
        reliability_mask            = squeeze(maxXCmat >= params.feature_params.minXCreliable);
%         reliability_mask            = squeeze(maxXCmat >= 0.9);
        azm_vec                     = params.global_params.azm_vec;
        num_azms                    = numel(azm_vec);
        peak_azm_vals               = zeros(num_glimpses,1);
        for zz = 1:num_glimpses
            [rowz,colz]             = find(glimpse_map==zz & reliability_mask==1);              % get row/column indices of the elements in the regarded glimpse
            if isempty(rowz)
                [rowz,colz]             = find(glimpse_map==zz);  
            end
            
            size_zz                 = length(rowz);                         % number of t-f-units units in the regarded glimpse
            acc_log_azs             = zeros(size_zz,num_azms);              % accumulate all log_like_az vectors beloging to this glimpse in one matrix
            for uu = 1:size_zz
                acc_log_azs(uu,:)   = feature_mat(:,rowz(uu),colz(uu));
            end
            mean_probs  = mean(acc_log_azs,1);                              % azimuth probabilities
            [~,Mc,~,~]           = findLocalPeaks(mean_probs.');            % get index of the peak azimuth for each glimpse
            if isempty(Mc)
                [~,Mc]           = max(mean_probs);                         % use the max value if no peak exists
            end
            peak_azm_vals(zz)  =  azm_vec(Mc);    
        end
        mean_glimpse_features{aa} = peak_azm_vals;                          % num_glimpses X 1

    elseif strcmp(feature_types{aa},'glimpse_size')==1                      % # of t-f-units in the glimpse
        glimpse_sizes = zeros(num_glimpses,1); 
        for zz = 1:num_glimpses
            [rowz,~]                 = find(glimpse_map==zz);               % get row and column indices of the elements in the regarded glimpse
            size_zz                  = length(rowz);                        % number of t-f-units in glimpse zz
            glimpse_sizes(zz)        = size_zz;
        end
        mean_glimpse_features{aa} = glimpse_sizes;                          % num_glimpses X 1
        
    elseif strcmp(feature_types{aa},'average_time')==1                      % average the indices of time frames to use temporal proximity        
        glimpse_avtime = zeros(num_glimpses,1);
        for zz = 1:num_glimpses
            [rowz,~]                  = find(glimpse_map==zz);              % get row and column indices of the elements in the regarded glimpse
            mean_rows                 = mean(rowz);             
            glimpse_avtime(zz)        = mean_rows*params.prep_params.hopsize_ms*0.001;
        end
        mean_glimpse_features{aa} = glimpse_avtime;
        
    elseif strcmp(feature_types{aa},'average_subband')==1                   % average the subband indices to be able to use spectral overlap/proximity
        glimpse_avsubband = zeros(num_glimpses,1); 
        for zz = 1:num_glimpses
            [~,colz]                = find(glimpse_map==zz);   
            mean_colz               = mean(colz);             
            glimpse_avsubband(zz)   = mean_colz; 
        end
        mean_glimpse_features{aa}   = glimpse_avsubband;
         
    elseif strcmp(feature_types{aa},'mean_subband_power')==1                % mean power in each frequency band (per glimpse)    
        feature_mat                 = 10*log10(feature_mats.summed_power_mat); % already summed over both ears
        glimpse_subband_power       = zeros(num_glimpses,params.prep_params.nFilts); % 1st col: number of the glimpse, 2nd: most likely f0, 3rd: height of the peak, 4th: # of t-f-units in the glimpse
        for zz = 1:num_glimpses
            [rowz,colz]             = find(glimpse_map==zz);    % get row and column indices of the elements in the regarded glimpse
            for bb = 1:params.prep_params.nFilts
                glimpse_subband_power(zz,bb) = mean(feature_mat(rowz(colz==bb),bb));  % mean power of t-f-units in subband bb in glimpse zz
            end
        end
        glimpse_subband_power(isnan(glimpse_subband_power))     = -25; % 
        glimpse_subband_power(find(glimpse_subband_power<-25))  = -25; % in case there are some really low energy filterbands, set them to the min value (could also be a lower value)
        mean_glimpse_features{aa}                               = glimpse_subband_power;
        
    elseif strcmp(feature_types{aa},'spectral_centroid')==1
        feature_mat                 = feature_mats.summed_power_mat;    % linear power, already summed over both ears
        spectral_centroids          = zeros(num_glimpses,1);
        for zz = 1:num_glimpses
            subband_powers              = zeros(params.prep_params.nFilts,1);
            [rowz,colz]                 = find(glimpse_map==zz);    % get row and column indices of the elements in the regarded glimpse
            for bb = 1:params.prep_params.nFilts
                subband_powers(bb,1) = mean(feature_mat(rowz(colz==bb),bb));  % mean power of t-f-units in subband bb in glimpse zz
            end
            subband_powers(isnan(subband_powers)) = 0; % if activity is 0 in the subband, the mean is NaN. Replace these entries with zeros
            spectral_centroids(zz) = sum(params.prep_params.cfs.*subband_powers.')/sum(subband_powers);
        end
        
        mean_glimpse_features{aa}                               = spectral_centroids;
        
    elseif strcmp(feature_types{aa},'fractional_subband_power')==1
       
        feature_mat                 = (feature_mats.summed_power_mat).^0.2; % already summed over both ears
        glimpse_subband_power       = zeros(num_glimpses,params.prep_params.nFilts); % 1st col: number of the glimpse, 2nd: most likely f0, 3rd: height of the peak, 4th: # of t-f-units in the glimpse
        for zz = 1:num_glimpses
            [rowz,colz]             = find(glimpse_map==zz);    % get row and column indices of the elements in the regarded glimpse
            for bb = 1:params.prep_params.nFilts              
                glimpse_subband_power(zz,bb) = mean(feature_mat(rowz(colz==bb),bb));  % mean power of t-f-units in subband bb in glimpse zz
            end
        end
        glimpse_subband_power(isnan(glimpse_subband_power))     = -25; % could also be a lower value
        glimpse_subband_power(find(glimpse_subband_power<-25))  = -25; % in case there are some really low energy filterbands, set them to the min value 
        mean_glimpse_features{aa}                               = glimpse_subband_power;
        
        
    elseif strcmp(feature_types{aa},'mean_power')==1              % averaged power in the glimpse        
        feature_mat                     = 10*log10(feature_mats.summed_power_mat); % this is already summed over both ears
        glimpse_pwrs = zeros(num_glimpses,1); % 1st col: number of the glimpse, 2nd: most likely f0, 3rd: height of the peak, 4th: # of t-f-units in the glimpse
        for zz = 1:num_glimpses
            [rowz,colz]                 = find(glimpse_map==zz);    % get row and column indices of the elements in the regarded glimpse
            size_zz                     = length(rowz);             % number of t-f-units in glimpse zz
            acc_pwr                     = zeros(size_zz, 1);  % OPTION at this point it might make more sense to first average within one time frame
            for uu = 1:size_zz                                  % accumulate all the periodicity vectors belonging to the regarded glimpse
                acc_pwr(uu)             = feature_mat(rowz(uu),colz(uu));
            end
            mean_pwr                    = mean(acc_pwr);    % sum over the first dimension to get index of most pronounced f0
            glimpse_pwrs(zz)            = mean_pwr;    % convert the index I to an actual f0
        end
        mean_glimpse_features{aa}       = glimpse_pwrs;
        
    elseif strcmp(feature_types{aa},'power_var')==1   % power variance (or here std) in each glimpse
        feature_mat                 = feature_mats.summed_power_mat; % already summed over both ears
        glimpse_pwr_std             = zeros(num_glimpses,1); % 1st col: number of the glimpse, 2nd: most likely f0, 3rd: height of the peak, 4th: # of t-f-units in the glimpse
        for zz = 1:num_glimpses
            [rowz,colz]                 = find(glimpse_map==zz);    % get row and column indices of the elements in the regarded glimpse
            size_zz                     = length(rowz);             % number of t-f-units in glimpse zz
            acc_pwr                     = zeros(size_zz, 1);  % OPTION at this point it might make more sense to first average within one time frame
            for uu = 1:size_zz                                  % accumulate all the periodicity vectors belonging to the regarded glimpse
                acc_pwr(uu)             = feature_mat(rowz(uu),colz(uu));
            end
            pwr_std                     = std(acc_pwr);    % sum over the first dimension to get index of most pronounced f0
            glimpse_pwr_std(zz)         = pwr_std;    % convert the index I to an actual f0
        end
        mean_glimpse_features{aa}   = glimpse_pwr_std;

     elseif strcmp(feature_types{aa},'peak_f0_val')==1  % average pitch (in Hz) in each glimpse
        feature_mat                 = feature_mats.complete_NAC_mat; 
        period_vec                  = params.feature_params.minLags : 1: params.feature_params.maxLags;
        f0_vec                      = (1./period_vec)*params.prep_params.fs; % fundamental frequencies in Hz
        num_f0s                     = numel(f0_vec);
        
        peak_f0_vals = zeros(num_glimpses,1); % 1st col: number of the glimpse, 2nd: most likely f0, 3rd: height of the peak, 4th: # of t-f-units in the glimpse
        for zz = 1:num_glimpses
            [rowz,colz]                 = find(glimpse_map==zz);    % get row and column indices of the elements in the regarded glimpse
            size_zz                     = length(rowz);             % number of t-f-units in glimpse zz
            acc_periodicity_vecs        = zeros(size_zz, num_f0s);  % OPTION at this point it might make more sense to first average within one time frame
            for uu = 1:size_zz                                  % accumulate all the periodicity vectors belonging to the regarded glimpse
                acc_periodicity_vecs(uu,:) = feature_mat(:,rowz(uu),colz(uu));
            end
            mean_periodicity_vec     = mean(acc_periodicity_vecs,1);    % average over the first dimension to get index of most pronounced f0
            
            [~,Mc,~,~]  = findLocalPeaks(mean_periodicity_vec.');            % get index of the peak period for each glimpse
            if isempty(Mc)
                [~,Mc]           = max(mean_periodicity_vec);                         % use the max value if no peak exists
            end
            
            peak_f0_vals(zz,:)      = f0_vec(Mc);
        end
        mean_glimpse_features{aa}       = peak_f0_vals; % num_glimpses X 1
        
    else
        error('undefined feature type')
    end
    
end


