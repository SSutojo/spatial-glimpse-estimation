function [ gmm_name ] = SVtypes2name( SV_types )
%SVVEC2GMM the SV vector goes in and the GMM name is created
%   Detailed explanation goes here
%   INPUT SV_types is a cell array. In the first column the feature type
%   and in the second column is the type of correlation coefficient
% 

num_SV_types = size(SV_types,1);

SV_cell = cell(num_SV_types,1); % slightly change the input cell SV_types in a way that feature type and correlation type are already in one string
for aa = 1:num_SV_types
    SV_cell{aa,1} = [SV_types{aa,1},'_',SV_types{aa,2}];    
end

name = [];

for idx = 1:num_SV_types
    
    test = 1;
    if strcmp(SV_cell{idx},'periodicity_degree_pearsons_corr' )==1
        name_component = 'PDPears';
    elseif strcmp(SV_cell{idx},'periodicity_degree_general_corr')==1
        name_component = 'PDGen';
    elseif strcmp(SV_cell{idx},'normalized_ac_pearsons_corr')==1
        name_component = 'NACPears';
    elseif strcmp(SV_cell{idx},'normalized_ac_general_corr')==1
        name_component = 'NACGen';        
    elseif strcmp(SV_cell{idx},'azimuth_probabilities_pearsons_corr')==1
        name_component = 'AZPears';
    elseif strcmp(SV_cell{idx},'log_azimuth_probabilities_pearsons_corr') == 1
        name_component = 'LOGAZPears';
    elseif strcmp(SV_cell{idx},'power_abs_diff') == 1
        name_component = 'PWRDelts';
    elseif strcmp(SV_cell{idx},'periodicity_degree_pearsons_corr_context' )==1
        name_component = 'PDPearsCon';
    elseif strcmp(SV_cell{idx},'periodicity_degree_general_corr_context')==1
        name_component = 'PDGenCon';
    elseif strcmp(SV_cell{idx},'normalized_ac_pearsons_corr_context')==1
        name_component = 'NACPearsCon';
    elseif strcmp(SV_cell{idx},'normalized_ac_general_corr_context')==1
        name_component = 'NACGenCon';
    elseif strcmp(SV_cell{idx},'normalized_ac_mean_max')==1
        name_component = 'NACMeanMax';    
    elseif strcmp(SV_cell{idx},'log_power_abs_diff')==1
        name_component = 'LOGPWRDelts';
%     elseif strcmp(SV_cell{idx},'log_azimuth_probabilities_pearsons_corr_context') == 1
%         name_component = 'LOGAZPearsCon';
%     elseif strcmp(SV_cell{idx},'power_abs_diff_context') == 1
%         name_component = 'PWRDeltsCon';
    elseif strcmp(SV_cell{idx},'power_abs_sum') == 1
        name_component = 'PWRSum';
    elseif strcmp(SV_cell{idx},'log_power_abs_sum') == 1
        name_component = 'LOGPWRSum';
    elseif strcmp(SV_cell{idx},'croot_power_abs_diff')==1
        name_component = 'CRTPWRDelts';
    elseif strcmp(SV_cell{idx},'croot_power_abs_sum') == 1
        name_component = 'CRTPWRSum';
    else
        error('unknown SV_type for pre-trained GMM')
    end
    
    name = [name,'',name_component];
    
    
end


gmm_name = name;



end

