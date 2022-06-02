function [featuresOUT] = normSingleFeatureDims(featuresIN,featureTypes,params)
%NORM_FEATURE_DIMS Normalize glimpse similarity features for global
%grouping (used in training and inference)
%
% the features are roughly normalized to a range from 0 to 1
% since this makes less sense for glimpse size and average time, these are
% normalized differently

% size featuresIN (nDatapoints X featureTypes)
% featureTypes [feature1A ; feature2A ; feature3A; ...; feature1B ; feature2B ; feature3B ; ...]
% feature1A is the first feature type in glimpse A


% featuresIN = featuresIN(:);


featuresOUT     = zeros(size(featuresIN));


if strcmp(featureTypes{1,1},'mean_power')==1           && strcmp(featureTypes{1,2},'mean')==1   % mean power
    featuresOUT = (featuresIN-min(featuresIN(:)))./(max(featuresIN(:))-min(featuresIN(:)));
    
elseif strcmp(featureTypes{1,1},'mean_power')==1       && strcmp(featureTypes{1,2},'abs_diff')==1   % mean power
    featuresOUT = (featuresIN-min(featuresIN(:)))./(max(featuresIN(:))-min(featuresIN(:)));
    
elseif strcmp(featureTypes{1,1},'mean_log_azm')==1     && strcmp(featureTypes{1,2},'peak_val')==1  % azm value
    featuresOUT = (featuresIN+90)./180;
    
elseif strcmp(featureTypes{1,1},'mean_xcorr')==1       && strcmp(featureTypes{1,2},'mean')==1   % mean power
    featuresOUT = (featuresIN-min(featuresIN(:)))./(max(featuresIN(:))-min(featuresIN(:)));
    
elseif strcmp(featureTypes{1,1},'mean_nac')==1         && strcmp(featureTypes{1,2},'peak_val')==1   % f0 value
    featuresOUT = (featuresIN-(params.prep_params.fs/params.feature_params.maxLags))./((params.prep_params.fs/params.feature_params.minLags)-(params.prep_params.fs/params.feature_params.maxLags));

elseif strcmp(featureTypes{1,1},'mean_nac')==1         && strcmp(featureTypes{1,2},'mean')==1  % pitch salience
    featuresOUT = (featuresIN(:)-min(featuresIN(:)))./(max(featuresIN(:))-min(featuresIN(:)));
    
elseif strcmp(featureTypes{1,1},'mean_nac')==1         && strcmp(featureTypes{1,2},'mean_max')==1  % pitch salience
    featuresOUT = (featuresIN-min(featuresIN(:)))./(max(featuresIN(:))-min(featuresIN(:)));
    
elseif strcmp(featureTypes{1,1},'mean_nac')==1         && strcmp(featureTypes{1,2},'std_max')==1  % difference of pitch salience
    featuresOUT = (featuresIN-min(featuresIN(:)))./(max(featuresIN(:))-min(featuresIN(:)));
    
elseif strcmp(featureTypes{1,1},'glimpse_size')==1     && strcmp(featureTypes{1,2},'direct')==1 % set in relation to number of T-F units in 0.01 sec
    ref_size    = (0.01*params.prep_params.nFilts)/(params.prep_params.hopsize_ms*0.001);
    currFeats   = featuresIN;
    currFeats(currFeats>ref_size) = ref_size;               % limit to a value of 1
    featuresOUT = (currFeats-1)./(ref_size);
    
elseif strcmp(featureTypes{1,1},'glimpse_size')==1     && strcmp(featureTypes{1,2},'max')==1 % set in relation to number of T-F units in 0.01 sec
    ref_size = (0.01*params.prep_params.nFilts)/(params.prep_params.hopsize_ms*0.001);
    currFeats = featuresIN;
    currFeats(currFeats>ref_size) = ref_size;               % limit to a value of 1
    featuresOUT(:,aa) = (currFeats-1)./(ref_size);
    
elseif strcmp(featureTypes{1,1},'glimpse_size')==1     && strcmp(featureTypes{1,2},'min')==1 % set in relation to number of T-F units in 0.01 sec
    ref_size = (0.01*params.prep_params.nFilts)/(params.prep_params.hopsize_ms*0.001);
    currFeats = featuresIN;
    currFeats(currFeats>ref_size) = ref_size;               % limit to a value of 1
    featuresOUT = (currFeats-1)./(ref_size);
    
elseif strcmp(featureTypes{1,1},'average_time')==1     && strcmp(featureTypes{1,2},'mean')==1 % set in relation to a common filelength of 4 sec
    featuresOUT = featuresIN./4;
    
elseif strcmp(featureTypes{1,1},'average_time')==1     && strcmp(featureTypes{1,2},'abs_diff')==1 % set in relation to a common filelength of 4 sec
    featuresOUT = featuresIN./4;
    
elseif strcmp(featureTypes{1,1},'average_subband')==1  && strcmp(featureTypes{1,2},'abs_diff')==1 % set in relation to a common filelength of 4 sec

    featuresOUT = 1-(featuresIN./params.prep_params.nFilts);
elseif strcmp(featureTypes{1,1},'power_var')==1  && strcmp(featureTypes{1,2},'abs_diff')==1 % set in relation to a common filelength of 4 sec
    featuresOUT = (featuresIN-min(featuresIN(:)))./(max(featuresIN(:))-min(featuresIN(:)));
elseif strcmp(featureTypes{1,1},'power_var')==1  && strcmp(featureTypes{1,2},'direct')==1 
    featuresOUT = (featuresIN-min(featuresIN(:)))./(max(featuresIN(:))-min(featuresIN(:)));

elseif strcmp(featureTypes{1,1},'mean_subband_power')==1  && strcmp(featureTypes{1,2},'direct')==1 % set in relation to a common filelength of 4 sec

    featuresOUT = (featuresIN-(-25))./(5-(-25));
elseif strcmp(featureTypes{1,1},'fractional_subband_power')==1  && strcmp(featureTypes{1,2},'direct')==1 % set in relation to a common filelength of 4 sec

    featuresOUT = (featuresIN-(-25))./(5-(-25));
    
elseif strcmp(featureTypes{1,1},'peak_azm_val')==1  && strcmp(featureTypes{1,2},'abs_diff')==1 % set in relation to a common filelength of 4 sec
        featuresOUT = 1-(featuresIN./180);      
elseif strcmp(featureTypes{1,1},'peak_azm_val')==1  && strcmp(featureTypes{1,2},'direct')==1 % set in relation to a common filelength of 4 sec
        featuresOUT = 1-(featuresIN./180); 
elseif strcmp(featureTypes{1,1},'reliable_peak_azm')==1  && strcmp(featureTypes{1,2},'direct')==1 % set in relation to a common filelength of 4 sec
        featuresOUT = 1-(featuresIN./180); 
               
        
elseif strcmp(featureTypes{1,1},'spectral_centroid')==1  && strcmp(featureTypes{1,2},'abs_diff')==1 % set in relation to a common filelength of 4 sec
 
        featuresOUT = 1-(featuresIN./params.prep_params.start_freqs(end));
elseif strcmp(featureTypes{1,1},'spectral_centroid')==1  && strcmp(featureTypes{1,2},'direct')==1 % set in relation to a common filelength of 4 sec

        featuresOUT = 1-(featuresIN./params.prep_params.start_freqs(end));
        
elseif strcmp(featureTypes{1,1},'peak_f0_val')==1         && strcmp(featureTypes{1,2},'abs_diff')==1   % f0 value
    featuresOUT = (featuresIN-(params.prep_params.fs/params.feature_params.maxLags))./((params.prep_params.fs/params.feature_params.minLags)-(params.prep_params.fs/params.feature_params.maxLags));
        

        
else
    featuresOUT = featuresIN;   % not all feature types need to be normalized as they are already in the 0-1 range
end


end

