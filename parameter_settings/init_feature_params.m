function [feature_params] = init_feature_params(p_features)
%INIT_FEATURE_PARAMS initialize the parameter struct for feature
%extraction, based on the independent parameters (set in
%def_feature_params), further (dependent) parameters are derived 

feature_params                  = p_features;


%% for binaural features
feature_params.azm_blocksize    = feature_params.fs *feature_params.azm_Lwin*0.001;

feature_params.maxDelay         = ceil(feature_params.maxITDsec*feature_params.fs); % same thing in samples
feature_params.lagsITD          = floor(-feature_params.maxITDsec*feature_params.fs):1:feature_params.maxDelay;

azm_classifier                  = load(['azmGMM_',num2str(feature_params.azm_Lwin),feature_params.azm_wintype,'.mat']);
feature_params.C                = azm_classifier.C;

%% for periodicity features
feature_params.pitch_blocksize  = feature_params.fs *feature_params.pitch_Lwin*0.001; % must be a whole-number multiple of the hopsize

feature_params.maxLags          = ceil(feature_params.fs/60);      % max lag for periodicity analysis in samples (z.B. max tau for 65Hz)
feature_params.minLags          = floor(feature_params.fs/400);    % min lag roughly 2.4 ms (420 Hz)


%% for summed_power features
feature_params.pow_blocksize    = feature_params.fs *feature_params.pow_Lwin*0.001;


end

