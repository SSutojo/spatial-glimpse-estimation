function [updated_params] = update_parameters_new(new_settings, prep_params, feature_params, local_params, global_params)
%UPDATE_PARAMETERS_NEW replace various fields in the params struct
%
%   INPUTs
%   new_settings: cell array, dimensions: 2 X number of updated parameters
%                 cells in first row contain the field names of parameters
%                 that are updated, cells in second row contain new values
%   prep_params, feature_params, local_params, global_params : old
%                 parameter structs
%
%   OUTPUTS
%   updated_params: the single parameter structs are joined and updated
%


params.prep_params = prep_params;
params.feature_params = feature_params;
params.local_params   = local_params;

if nargin > 4
    params.global_params = global_params;
end




% get number of parameters to update
num_updates = size(new_settings,2);
for ii = 1:num_updates
    if isempty(new_settings{2,ii}) == 1
        h = 1;
        continue
    elseif isempty(new_settings{2,ii}) ~= 1
        hh = 1;
        % get the position of the field in the struct and the fieldname
        pathlength = size(new_settings{1,ii},2);
        if pathlength == 1 % find a better way to do this someday
            params.(new_settings{1,ii}) = new_settings{2,ii}; % replace by new value
        elseif pathlength == 2
            params.(new_settings{1,ii}{1,1}).(new_settings{1,ii}{1,2}) = new_settings{2,ii};
        elseif pathlength == 3
            params.(new_settings{1,ii}{1,1}).(new_settings{1,ii}{1,2}).(new_settings{1,ii}{1,3}) = new_settings{2,ii};
        elseif pathlength == 4
            params.(new_settings{1,ii}{1,1}).(new_settings{1,ii}{1,2}).(new_settings{1,ii}{1,3}).(new_settings{1,ii}{1,4}) = new_settings{2,ii};
        elseif pathlength == 5
            params.(new_settings{1,ii}{1,1}).(new_settings{1,ii}{1,2}).(new_settings{1,ii}{1,3}).(new_settings{1,ii}{1,4}).(new_settings{1,ii}{1,5}) = new_settings{2,ii};
        elseif pathlength == 6
            params.(new_settings{1,ii}{1,1}).(new_settings{1,ii}{1,2}).(new_settings{1,ii}{1,3}).(new_settings{1,ii}{1,4}).(new_settings{1,ii}{1,5}).(new_settings{1,ii}{1,6}) = new_settings{2,ii};
        else
            error('check updated fieldname')  % kann man auch einfach erweitern, bisher hab ich nicht mehr als 6 unterstructs
        end
    end
    
    
end



% the new settings must contain the name of the field and the new value

% re initialize the parameter structs
updated_params.prep_params     = init_prep_params(params.prep_params);
updated_params.feature_params  = init_feature_params(params.feature_params);
updated_params.local_params    = init_local_params(params.local_params);

if isfield(params,'global_params') == 1
    updated_params.global_params = init_global_params(params.global_params);
end





end

