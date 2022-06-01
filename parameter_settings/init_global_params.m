function [global_params] = init_global_params(globpars)
%INIT_GLOBAL_PARAMS initialize the parameter struct for global grouping, based
%on the independent parameters (set in def_global_params), further (dependent) parameters are derived 

global_params = globpars;
[global_params.feature_types]       = glob_control_order(global_params.feature_types);  % make sure features are in the correct order

% Select optimizer
global_params.model_confs.network.optimizer = global_params.model_confs.network.optimizer.(upper(global_params.model_confs.network.optimizer.algorithm));

end

