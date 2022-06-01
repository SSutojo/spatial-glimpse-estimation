function [def_global] = def_global_params
% define independent parameters for global grouping (grouping glimpses
% over a larger time span to obtain an IBM estimate)


def_global.est_num_sources          = 4;   % estimated number of sources (if true number of sources is known, this field can be replaced later)

def_global.group_by                 = 'absolute_glimpse_features'; % 'relative_glimpse_features','stacked_glimpse_features'

def_global.grouping_method          = 'histogram1D';             
def_global.feature_types            = {'reliable_peak_azm','direct'}; % feature types 


def_global.hist_peaks               = 'per_time_frame';     % 'per_t_f_unit', 'per_time_frame', 'per_glimpse'

def_global.SNRcrit                  = 10;                                                
                                                                                        
%% parameters for each method
def_global.min_glimpse_size     = 1; % minimum glimpse size in T-F units, option to consider only larger glimpses in the source identification
                                     % if all glimpses should be considered, set this to 1
%********** periodicity (pitch)
def_global.low_f0_emphasis      = true; % if this is 'true', the higher frequencies are less weighted when detecting prominent f0s
def_global.lowpass_fc           = 400;   % cutoff frequency for the lowpass to emphasise lower freqs
def_global.pitchsal_weight      = true;  % use pitch salience as weight
def_global.smooth_featrates     = true;  % apply gaussian smoothing to the f0 axis before determining the peaks/max values

def_global.pitch_tolerance_radius    = 30;  % within this radius around a peak, 30 is good for periodicity only ---vorher 30
                                                % t-f- units are assigned to the source (here that's 20 degree azimuth)
                                                
%********** azimuth (position)
def_global.azm_tolerance_radius      = 5;  % vorher 5
def_global.xcorr_weigth              = false;
def_global.azm_vec                   = -90:5:90; 



%% Parameters for DNN training / glimpse tracker training
def_global.model_confs.model_type              = 'dnn';
def_global.model_confs.label_type              = 'binary';    % *v* soft or binary
% These parameters do not affect the feature space that is stored
def_global.model_confs.post.label = 'post-processing parameters';
def_global.model_confs.post.trainRatio = 1;   % 1 if all

% Cell array of features used for training
def_global.model_confs.post.strFeatures = {'fbe'};
def_global.model_confs.post.idxFeatures = {[]};
def_global.model_confs.post.compress = {'cuberoot'};

% Feature normalization
def_global.model_confs.post.normalize.global.normMethod        = 'none'; % *v*    can be 'none' 'mean' 'var' 'meanvar' 'max' 'range' 'heq' -- vorher 'none'
def_global.model_confs.post.normalize.local.normMethod         = 'heq'; % *v*    can be 'heq', 'heq_overchans'(same norm for all freq channels),'meanvar_overchans'

% Context expansion
def_global.model_confs.post.context.context                    = [0 0]; % *v*
def_global.model_confs.post.context.bRemoveDC                  = false;
def_global.model_confs.post.context.bDCT                       = false;
def_global.model_confs.post.context.orderDCT                   = 5;
def_global.model_confs.post.context.bAppend                    = false;

%% NETWORK PARAMETERS
def_global.model_confs.network.label = 'neural network configuration';

% Network architecture
def_global.model_confs.network.architecture.nLayersHidden      = 3;   % *v*  3
% opt.network.architecture.sizeLayersHidden   = size(global_params.grouping_features,1); % *v* 10    mit 2 multiplizieren
def_global.model_confs.network.architecture.sizeLayersHidden   = 4; % *v* 10    mit 2 multiplizieren
def_global.model_confs.network.architecture.activationHidden   = 'relu';  % relu
def_global.model_confs.network.architecture.activationOut      = 'sigmoid'; % sigmoid
def_global.model_confs.network.architecture.nEnsemble          = 1;

% Network regularization
def_global.model_confs.network.regularization.lambdaL1Norm     = 0;
def_global.model_confs.network.regularization.lambdaL2Norm     = 0;
def_global.model_confs.network.regularization.maxL2Norm        = 0;  % *v* 0
def_global.model_confs.network.regularization.dropoutInput     = 0;
def_global.model_confs.network.regularization.dropoutHidden    = 0.2; % *v* vorher 0.2
def_global.model_confs.network.regularization.bInvertedDropout = true;

def_global.model_confs.network.regularization.bDropoutLast     = false;
def_global.model_confs.network.regularization.bBatchNorm       = true; % *v*
def_global.model_confs.network.regularization.bBNBeforeActivation = true; % *v*
def_global.model_confs.network.regularization.alphaBN          = 0.99;
def_global.model_confs.network.regularization.epsilonBN        = 0.001;

% Network training
def_global.model_confs.network.training.bSaveMemory            = false;
def_global.model_confs.network.training.bSpeedUp               = false;
def_global.model_confs.network.training.initWeights            = 'he_uniform'; % vorher glorot

def_global.model_confs.network.training.initBias               = 'weights'; % 'weights' = same as weights, or a scalar
def_global.model_confs.network.training.nEpochs                = 300; % vorher 300
def_global.model_confs.network.training.earlyStopping          = 25; % vorher 20
def_global.model_confs.network.training.batchSize              = 256; % vorher 128

def_global.model_confs.network.training.metric                 = 'loss'; % either 'loss' or 'performance'
def_global.model_confs.network.training.lossFunction           = 'mse'; % 'mse' 'mae' 'msle' 'ce' ...
def_global.model_confs.network.training.performanceFunction    = 'hfa'; % 'hfa'

def_global.model_confs.network.training.learningRate           = 0.001; % before 0.001
def_global.model_confs.network.training.learningRateSchedule   = 'none'; % 'none' 'step' '1/t' 'exp'
def_global.model_confs.network.training.learningRateDecay      = 0.9;
def_global.model_confs.network.training.learningRateDecayStep  = 5;

def_global.model_confs.network.training.labelShift             = false;  % if true, the labels are multiplied with 0.7 and 0.15 is added

% Diagnostics during training
def_global.model_confs.network.diagnostics.bPlot               = true;
def_global.model_confs.network.diagnostics.bVerbose            = true;
def_global.model_confs.network.diagnostics.plotStep            = 5;
def_global.model_confs.network.diagnostics.bWeightHist         = true;
def_global.model_confs.network.diagnostics.bWeightRatio        = true;
def_global.model_confs.network.diagnostics.rangeHist           = [-1 1];
def_global.model_confs.network.diagnostics.nBinsHist           = 100;

% Network optimization
def_global.model_confs.network.optimizer.algorithm             = 'adam2';  %model_confs.opt.network.optimizer.algorithm = P.optimizer;

% Optimizer-specific settings
def_global.model_confs.network.optimizer.SGD.label             = 'sgd';
def_global.model_confs.network.optimizer.SGD                   = [];

def_global.model_confs.network.optimizer.SGD_MOM.label         = 'sgd_mom';
def_global.model_confs.network.optimizer.SGD_MOM.switchEpoch   = 5;
def_global.model_confs.network.optimizer.SGD_MOM.decayInitial  = 0.5;
def_global.model_confs.network.optimizer.SGD_MOM.decayFinal    = 0.9;

def_global.model_confs.network.optimizer.NAG.label             = 'nag';
def_global.model_confs.network.optimizer.NAG.switchEpoch       = 5;
def_global.model_confs.network.optimizer.NAG.decayInitial      = 0.5;
def_global.model_confs.network.optimizer.NAG.decayFinal        = 0.9;

def_global.model_confs.network.optimizer.ADA_DELTA.label       = 'ada_delta';
def_global.model_confs.network.optimizer.ADA_DELTA.decay       = 0.9;
def_global.model_confs.network.optimizer.ADA_DELTA.epsilon     = 1E-8;

def_global.model_confs.network.optimizer.ADA_MAX.label         = 'ada_max';
def_global.model_confs.network.optimizer.ADA_MAX.decay1        = 0.9;
def_global.model_confs.network.optimizer.ADA_MAX.decay2        = 0.999;
def_global.model_confs.network.optimizer.ADA_MAX.epsilon       = eps;

def_global.model_confs.network.optimizer.RMS_PROP.label        = 'rms_prop';
def_global.model_confs.network.optimizer.RMS_PROP.decay        = 0.9;
def_global.model_confs.network.optimizer.RMS_PROP.epsilon      = 1E-8;

def_global.model_confs.network.optimizer.ADAM.label            = 'adam';
def_global.model_confs.network.optimizer.ADAM.decay1           = 0.9;
def_global.model_confs.network.optimizer.ADAM.decay2           = 0.999;
def_global.model_confs.network.optimizer.ADAM.epsilon          = 1E-8;

def_global.model_confs.network.optimizer.ADAM2.label           = 'adam2';
def_global.model_confs.network.optimizer.ADAM2.decay1          = 0.9;
def_global.model_confs.network.optimizer.ADAM2.decay2          = 0.999;
def_global.model_confs.network.optimizer.ADAM2.epsilon         = 1E-8;



% training data
def_global.model_confs.training_data.MixTypes      = '2sp,3sp,4sp';
def_global.model_confs.training_data.NoiseType     = 'Low,Medium';
def_global.model_confs.training_data.Size          = 'test';
def_global.model_confs.training_data.Balancing     = 'none';








