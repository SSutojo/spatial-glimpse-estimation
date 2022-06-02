function [out_cell] = glob_control_order(in_cell)
%GLOB_CONTROL_ORDER ensure that the feature types for global grouping are
%always in the same order

num_feature_types   = size(in_cell,1);
rank_cell           = cell(num_feature_types,3);
rank_cell(:,2:3)    = in_cell;

for aa = 1:num_feature_types
    if strcmp(rank_cell{aa,2},'peak_azm_val')==1    && strcmp(rank_cell{aa,3},'direct')==1
        rank_cell{aa,1}         = 1; 
    elseif strcmp(rank_cell{aa,2},'peak_azm_val')==1    && strcmp(rank_cell{aa,3},'abs_diff')==1
        rank_cell{aa,1}         = 2; 
    elseif strcmp(rank_cell{aa,2},'mean_subband_power')==1    && strcmp(rank_cell{aa,3},'direct')==1
        rank_cell{aa,1}         = 3; 
    elseif strcmp(rank_cell{aa,2},'spectral_centroid')==1    && strcmp(rank_cell{aa,3},'abs_diff')==1
        rank_cell{aa,1}         = 4; 
    elseif strcmp(rank_cell{aa,2},'spectral_centroid')==1    && strcmp(rank_cell{aa,3},'direct')==1
        rank_cell{aa,1}         = 5; 
    elseif strcmp(rank_cell{aa,2},'power_var')==1    && strcmp(rank_cell{aa,3},'direct')==1
        rank_cell{aa,1}         = 6; 
    elseif strcmp(rank_cell{aa,2},'power_var')==1    && strcmp(rank_cell{aa,3},'abs_diff')==1
        rank_cell{aa,1}         = 7;      
    elseif strcmp(rank_cell{aa,2},'mean_power')==1    && strcmp(rank_cell{aa,3},'abs_diff')==1
        rank_cell{aa,1}         = 8;  
    elseif strcmp(rank_cell{aa,2},'peak_f0_val')==1    && strcmp(rank_cell{aa,3},'abs_diff')==1
        rank_cell{aa,1}         = 9;   
    elseif strcmp(rank_cell{aa,2},'reliable_peak_azm')==1    && strcmp(rank_cell{aa,3},'direct')==1
        rank_cell{aa,1}         = 10;   
    else
        error('unknown SV type')
    end     
end 
                                
ranked_cell         = sortrows(rank_cell,1);
out_cell            = ranked_cell(:,2:3);

end

