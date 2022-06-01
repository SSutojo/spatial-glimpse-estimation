function [def_features] = def_feature_params
% define the independent parameters for the extraction of auditory features

def_features.fs           = 16000;      % sampling frequency [Hz]
def_features.hopsize_ms   = 10;         % hopsize between time frames [ms]. frame lengths must be multiples of the hopsize

% length of time frames
def_features.azm_Lwin   = 20;           % length of time frames for azimuth estimation [ms]
def_features.pitch_Lwin = 30;           % length of time frames for pitch estimation [ms]
def_features.pow_Lwin   = 30;           % length of time frames for power feature estimation [ms]

% for binaural features
def_features.azm_wintype      = 'rect'; % window type for azimuth estimation (other options 'rect', 'hann' or 'shann')
def_features.maxITDsec        = 1.1e-3; % max ITD that's tested [s]
def_features.minXCreliable    = 0.3;    
def_features.ITDsec           = 1.1e-3; % min ITD that's tested [s]
def_features.TDOAsec          = 1e-4;
def_features.removeEndpoints  = true;
def_features.HaircellModel    = 'none'; % only used, if the "roman" haircell model is used for azimuth classifier


% for pitch features
def_features.pitch_wintype    = 'rect'; % window type for pitch estimation (other options 'rect', 'hann' or 'shann')
def_features.pitch_scale      = 'coeff';          % options: 'biased','unbiased','none','coeff'
def_features.model            = 'haircell_yes';   % option to use haircell model for pitch estimation

% for power features
def_features.pow_wintype        = 'rect'; % window type for pitch estimation (other options 'rect', 'hann' or 'shann')
