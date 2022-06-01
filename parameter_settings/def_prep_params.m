function [def_prep] = def_prep_params
% define independent parameters for pre processing and stimulus generation

def_prep.fs                  = 16000;       % sampling frequency [Hz]
def_prep.nChans              = 2;           % nr. of input channels 

%% stimulus generation
def_prep.HRIR_type           = 'KEMAR';                     % HRIR database 
def_prep.diffuseDirs         = [-175:5:180 ; zeros(1,72)];  % [-175:5:180 ; zeros(1,72)]
def_prep.noise               = 'diffuse_pink';              % background noise (options: 'white','diffuse','icra5','diffuse_pink')

def_prep.CARL                = false;                       % true if running on CARL cluster

%% temp segmentation (hopsize must be the same for every feature to maintain the same number of frames across all feature domains)
def_prep.hopsize_ms          = 10;                          % hopsize between time frames [ms]. Note: 
def_prep.Lwin_min            = def_prep.hopsize_ms *2;      % frame length for ground truth extraction [ms]
def_prep.window_type         = 'shan';                      % window type for ground truth extraction 'rect','shan'


%% filterbank (decomposition and resynthesis)
def_prep.f_low               = 200;   % lowest center frequency [Hz]
def_prep.f_high              = 7000;  % highest center frequency [Hz]
def_prep.filters_per_ERB     = 1.25;  % filter bandwidth [1/ERB]

def_prep.fb_type             = 'am_t';  % filterbank type (options: 'may','am_t')

%% haircell model
def_prep.haircell_type       = 'halfwave_smooth';  % type of haircell model
% required if 'halfwave_smooth' is used
def_prep.smooth_order       = 4;                   % order of the lowpass butterworth filter used to smooth the envelope
def_prep.smooth_freq        = 1000;                % cutoff frequency of the smoothing filter



