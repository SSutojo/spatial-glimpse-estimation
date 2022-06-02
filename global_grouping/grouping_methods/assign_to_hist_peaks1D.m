function [overall_feature_peaks, feature_rates, glimpse_ID_vec] = assign_to_hist_peaks1D(feature_mat,glimpse_features,options)
%ASSIGN_TO_HIST_PEAKS get the distribution/rates of feature values across the entire T-F plane to
%determine the most promiment feature values (e.g. the most common spatial
%location). Choose the most promiment feature values as belonging to the
%estimated sources (e.g. source locations) and assign the individual T-F units
%to these estimated sources
%
% IN 
%   feature_mat     - features per time-fraquency unit
%                   (dims: length of feature vec X num. time frames X num. frequency bands)
%   glimpse features - features per glimpse
%                   (dims: num glimpses X length of feature vec)
%   options         
%     .feature_vals - actual feature values that the feature can take, i.e. that correspond to the elements in the feature vector 
%                     (e.g. the azimuth values from -90:5:90 degree correspond to the 37 elements in the feature azm feature vector)
%     .hist_peaks   - defines if the histogram should be based on the
%                     feature peak in either individual T-F unit,
%                     time-frame, or glimpse
%     .est_num_sources - estimated number of sources
%     .tolerance_radius -  minimum distance between two peaks and radius
%                     around one peak in which every glimpse is assigned
%                     to this peak
%                   

num_glimpses                    = size(glimpse_features,1);              % get the number of separate clusters
num_win                         = size(feature_mat,2);
num_bands                       = size(feature_mat,3);
num_feature_vals                = numel(options.feature_vals);


%% source detection
% calculate the most pronounced azimuths across entire time-frequency plane,
% 3 possible methods: either sum across all t-f-units (method 1), or first determine the peak
% feature within one time frame to create a histogram (method 2), or use the peak feature 
% within each glimpse to create the histogram (method 3) - then, pick the
% highest peaks in this histogram

if strcmp(options.hist_peaks,'per_t_f_unit')==1
    max_feature_inds                = zeros(num_win, num_bands); 
    for aa = 1:num_win                                 
        for bb = 1:num_bands
            [~,M,~,~]           = findLocalPeaks(feature_mat(:,aa,bb)); % get the peak azimuth for each t-f-unit
            if isempty(M)
                [~,M]           = max(feature_mat(:,aa,bb));            % use the max value if no peak exists
            end
            max_feature_inds(aa,bb) = M;
        end
    end
    max_feature_inds                = max_feature_inds(:);
    feature_rates                   = zeros(num_feature_vals,2);
    for cc = 1:num_feature_vals
        feature_rates(cc,1)         = options.feature_vals(cc);               % determine the most frequent azimuths over entire t-f-plane
        feature_rates(cc,2)         = sum(max_feature_inds==cc);              % count the number of T-F units with the same index (same max feat)
    end
    
elseif strcmp(options.hist_peaks,'per_time_frame')==1                  % method 2
    subband_averaged            = sum(feature_mat,3);                   % sum log azm probs across all subbands in one time frame
    max_feature_inds            = zeros(num_win,1);
    for xx = 1:num_win
        [~,M,~,~]               = findLocalPeaks(subband_averaged(:,xx)); % or [~, max_inds]= max(subband_averaged,[],1);% get max per time_frame
        if isempty(M)
            [~,M]               = max(subband_averaged(:,xx));
        end
        max_feature_inds(xx)            = M;
    end
    feature_rates                   = zeros(num_feature_vals,2);
    for cc = 1:num_feature_vals                                              % count how many times a max_azm appears (could also use histcount or so)                         % 1st col: 
        feature_rates(cc,1)         = options.feature_vals(cc); 
        feature_rates(cc,2)         = sum(max_feature_inds == cc);                    % 2nd col: num. of time frames in which this azimuth has highest probability
    end
    
elseif strcmp(options.hist_peaks,'per_glimpse')==1
    feature_rates   = zeros(num_feature_vals,2);
    for cc = 1:num_feature_vals
        feature_rates(cc,1) = options.feature_vals(cc);
        feature_rates(cc,2) = sum(glimpse_features == options.feature_vals(cc));
    end

else
    error('define the method for determining the overall feature peaks')
end

%% choose peaks in the histogram (the signify the feature vals that occur most often)
[~ , sorted_inds]               = sort(feature_rates(:,2),'descend');       % get feature values with the highest hist_counts
overall_feature_peaks           = repmat((options.est_num_sources-1),1);   % 
overall_feature_peaks(1)        = options.feature_vals(sorted_inds(1));

% avoid that the overall peaks lie too close to each other (closer than one tolerance radius)
if options.est_num_sources >= 3
    counter                         = 1;
    for bb = 2:(options.est_num_sources-1)          % find peak features for est_num_sources-1 and make the last source a trash class
        previous_peak_features      = overall_feature_peaks;
        overall_feature_peaks(bb)   = options.feature_vals(sorted_inds(1+counter));
        while any(abs(overall_feature_peaks(bb)-previous_peak_features) < options.tolerance_radius) % in case the new peak is too close to another peak ...
            counter                 = counter+1;                            % ... just take the next azm with a high prob
            overall_feature_peaks(bb) = options.feature_vals(sorted_inds(1+counter));
        end
        counter                 = counter + 1;
    end
end

%% assign each glimpse to one of the histogram peaks, if it's not close to any peak, assign it to the trash class
glimpse_ID_vec              = zeros(num_glimpses,1);                
for pp = 1:num_glimpses                                             % for each glimpse: check if it is near an overall_feature_peak
    for qq = 1:(options.est_num_sources-1)                          % if so: assign the glimpse to the peaks' location
        curr_peak_angle     = overall_feature_peaks(qq);
        if (curr_peak_angle-options.tolerance_radius) <= glimpse_features(pp) && glimpse_features(pp) <= curr_peak_angle+options.tolerance_radius
            glimpse_ID_vec(pp)      = qq;                           % assign the glimpse an ID that now corresponds to this location
        else
            continue                                                % glimpse_ID_vec(pp) remains zero and is assigned to the 'leftovers' class
        end
    end
    if glimpse_ID_vec(pp)==0                                        % if this glimpse was not yet assigned any of the histogram peaks
        glimpse_ID_vec(pp)      = options.est_num_sources;          % ... it is assigned to the 'trash class'(leftovers)
    end
end

end

