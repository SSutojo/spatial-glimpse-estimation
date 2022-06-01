function [local_params] = init_local_params(Plocal)
%INIT_LOCAL_PARAMS initialize the parameter struct for local grouping, based
%on the independent parameters (set in def_prep_params), further (dependent) parameters are derived 


local_params = Plocal;

% control SV order
local_params.SV_types = SVs_control_order(local_params.SV_types);

% Select optimizer
local_params.model_confs.network.optimizer = local_params.model_confs.network.optimizer.(upper(local_params.model_confs.network.optimizer.algorithm));



end

