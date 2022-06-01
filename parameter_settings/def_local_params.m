function [def_local] = def_local_params
% define the independent parameters for local clustering (glimpse formation
% and estimation of feature contrasts)

% feature types for local clustering (relative features)
def_local.SV_types   = {'normalized_ac','mean_max';'normalized_ac','pearsons_corr';'log_power', 'abs_diff';'log_power', 'abs_sum';'azimuth_probabilities','pearsons_corr';'log_azimuth_probabilities','pearsons_corr'}; % 

% method for glimpse formation
def_local.contour_method                    = 'DilationErosion';
def_local.regiongrow.contrast_thresh        = 0.08; % used for regiongrow(1st value only), boundary_exp(1st value only), soft_boundary_exp(all values)

def_local.gb_superpixels.tau                = 0.13; % between 0 and 1. higher factor yields larger clusters

def_local.slic_superpixels.slic_compactness = 1;  % 10 for more compact, i.e. more square shaped
def_local.slic_superpixels.slic_stepsize    = 4;  % stepsize with which superpixel centers are distributed across image

def_local.rep_bound.contrast_threshs        = [0.05, 0.15, 0.35];
def_local.rep_bound.Nexpansions             = 1; % start iteration must be higher (more conservative) than end iteration

def_local.soft_bound.start_end_iteration    = [3 0];  % number of expansions per criterion in soft_boundary_expansion
def_local.soft_bound.contrast_thresh        = 0.3;  % decide of thresh1 ist start should be higher than thresh2 (end)


%% settings for the model to derive local contrasts
def_local.model_confs.model_type            = 'dnn';        % 'gmm' or 'dnn' or 'none'
def_local.model_confs.label_type            = 'soft';       % options: 'soft','binary'
def_local.model_confs.model_range           = 'broadband_sep'; %  'subband_sep' model per filterband, 'broadband_sep' model for all bands (across separate from within band transitions
                                                %  or 'broadband_all'(combines across+within frequency band transitions, most comprehensive model)

%% FEATURE POST-PROCESSING

% These parameters do not affect the feature space that is stored
def_local.model_confs.post.label                              = 'post-processing parameters';
def_local.model_confs.post.trainRatio                         = 1;   % 1 if all

% Cell array of features used for training
def_local.model_confs.post.strFeatures                        = {'fbe'};
def_local.model_confs.post.idxFeatures                        = {[]};
def_local.model_confs.post.compress                           = {'cuberoot'};

% Feature normalization
def_local.model_confs.post.normalize.global.normMethod        = 'none'; % option: 'none' 'mean' 'var' 'meanvar' 'max' 'range' 'heq'
def_local.model_confs.post.normalize.local.normMethod         = 'heq'; % options: 'heq', 'heq_overchans'(same norm for all freq channels),'meanvar_overchans'

% Context expansion
def_local.model_confs.post.context.context                    = [2 2]; % *v*
def_local.model_confs.post.context.bRemoveDC                  = false;
def_local.model_confs.post.context.bDCT                       = false;
def_local.model_confs.post.context.orderDCT                   = 5;
def_local.model_confs.post.context.bAppend                    = false;

%% NETWORK PARAMETERS
def_local.model_confs.network.label = 'neural network configuration';

% Network architecture
def_local.model_confs.network.architecture.nLayersHidden      = 2;   
def_local.model_confs.network.architecture.sizeLayersHidden   = 960; 
def_local.model_confs.network.architecture.activationHidden   = 'relu';
def_local.model_confs.network.architecture.activationOut      = 'sigmoid';
def_local.model_confs.network.architecture.nEnsemble          = 1;

% Network regularization
def_local.model_confs.network.regularization.lambdaL1Norm     = 0;
def_local.model_confs.network.regularization.lambdaL2Norm     = 0;
def_local.model_confs.network.regularization.maxL2Norm        = 0;  
def_local.model_confs.network.regularization.dropoutInput     = 0;
def_local.model_confs.network.regularization.dropoutHidden    = 0.2;
def_local.model_confs.network.regularization.bInvertedDropout = true;

def_local.model_confs.network.regularization.bDropoutLast     = false;
def_local.model_confs.network.regularization.bBatchNorm       = true; 
def_local.model_confs.network.regularization.bBNBeforeActivation = true; 
def_local.model_confs.network.regularization.alphaBN          = 0.99;
def_local.model_confs.network.regularization.epsilonBN        = 0.001;

% Network training
def_local.model_confs.network.training.bSaveMemory            = false;
def_local.model_confs.network.training.bSpeedUp               = false;
def_local.model_confs.network.training.initWeights            = 'he_uniform'; % vorher glorot

def_local.model_confs.network.training.initBias               = 'weights'; % 'weights' = same as weights, or a scalar
def_local.model_confs.network.training.nEpochs                = 300; % vorher 300
def_local.model_confs.network.training.earlyStopping          = 20;
def_local.model_confs.network.training.batchSize              = 128;

def_local.model_confs.network.training.metric                 = 'loss'; % 'loss','performance'
def_local.model_confs.network.training.lossFunction           = 'pow'; % 'mse' 'mae' 'msle' 'ce' ...
def_local.model_confs.network.training.performanceFunction    = 'hfa';

def_local.model_confs.network.training.learningRate           = 0.001;
def_local.model_confs.network.training.learningRateSchedule   = 'none'; % 'none' 'step' '1/t' 'exp'
def_local.model_confs.network.training.learningRateDecay      = 0.9;
def_local.model_confs.network.training.learningRateDecayStep  = 5;

def_local.model_confs.network.training.labelShift             = false;  
% Diagnostics during training
def_local.model_confs.network.diagnostics.bPlot               = false;
def_local.model_confs.network.diagnostics.bVerbose            = true;
def_local.model_confs.network.diagnostics.plotStep            = 5;
def_local.model_confs.network.diagnostics.bWeightHist         = true;
def_local.model_confs.network.diagnostics.bWeightRatio        = true;
def_local.model_confs.network.diagnostics.rangeHist           = [-1 1];
def_local.model_confs.network.diagnostics.nBinsHist           = 100;

% Network optimization
def_local.model_confs.network.optimizer.algorithm             = 'adam2';  %model_confs.opt.network.optimizer.algorithm = P.optimizer;

% Optimizer-specific settings
def_local.model_confs.network.optimizer.SGD.label             = 'sgd';
def_local.model_confs.network.optimizer.SGD                   = [];

def_local.model_confs.network.optimizer.SGD_MOM.label         = 'sgd_mom';
def_local.model_confs.network.optimizer.SGD_MOM.switchEpoch   = 5;
def_local.model_confs.network.optimizer.SGD_MOM.decayInitial  = 0.5;
def_local.model_confs.network.optimizer.SGD_MOM.decayFinal    = 0.9;

def_local.model_confs.network.optimizer.NAG.label             = 'nag';
def_local.model_confs.network.optimizer.NAG.switchEpoch       = 5;
def_local.model_confs.network.optimizer.NAG.decayInitial      = 0.5;
def_local.model_confs.network.optimizer.NAG.decayFinal        = 0.9;

def_local.model_confs.network.optimizer.ADA_DELTA.label       = 'ada_delta';
def_local.model_confs.network.optimizer.ADA_DELTA.decay       = 0.9;
def_local.model_confs.network.optimizer.ADA_DELTA.epsilon     = 1E-8;

def_local.model_confs.network.optimizer.ADA_MAX.label         = 'ada_max';
def_local.model_confs.network.optimizer.ADA_MAX.decay1        = 0.9;
def_local.model_confs.network.optimizer.ADA_MAX.decay2        = 0.999;
def_local.model_confs.network.optimizer.ADA_MAX.epsilon       = eps;

def_local.model_confs.network.optimizer.RMS_PROP.label        = 'rms_prop';
def_local.model_confs.network.optimizer.RMS_PROP.decay        = 0.9;
def_local.model_confs.network.optimizer.RMS_PROP.epsilon      = 1E-8;

def_local.model_confs.network.optimizer.ADAM.label            = 'adam';
def_local.model_confs.network.optimizer.ADAM.decay1           = 0.9;
def_local.model_confs.network.optimizer.ADAM.decay2           = 0.999;
def_local.model_confs.network.optimizer.ADAM.epsilon          = 1E-8;

def_local.model_confs.network.optimizer.ADAM2.label           = 'adam2';
def_local.model_confs.network.optimizer.ADAM2.decay1          = 0.9;
def_local.model_confs.network.optimizer.ADAM2.decay2          = 0.999;
def_local.model_confs.network.optimizer.ADAM2.epsilon         = 1E-8;



% training data
def_local.model_confs.training_data.MixTypes      = '2sp,3sp,4sp';
def_local.model_confs.training_data.NoiseType     = 'Low,Same,High';
def_local.model_confs.training_data.Size          = 'large';
def_local.model_confs.training_data.Balancing     = 'none';
def_local.model_confs.training_data.Rooms         = {'anechoic'};
% def_local.model_confs.training_data.Rooms         = {'anechoic','290ms','480ms','690ms'};



