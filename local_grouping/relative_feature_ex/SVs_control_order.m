function [ out_cell ] = SVs_control_order( in_cell )
%SVS_CONTROL_ORDER for consistency between the extracted features
%and the pre-trained GMMs, the order of the different Similarity Value (SV)
%types must stay the same throughout training and testing. To ensure the
%correct order,this function ranks the SVs types
% 
%   INPUTS: in_cell     - N X 2 cell with N being the number of SV types
%                         every row is an SV type,meaning a combination of feature
%                         type (1st column) and correlation type (2nd column)
%                         (in_cell is usually the 'local clustering mode')
%   OUTPUTS: out_cell   - cell of the same size as in_cell but the rows are
%                         swapped to match the desired order


num_SV_types        = size(in_cell,1);
rank_cell           = cell(num_SV_types,3);
rank_cell(:,2:3)    = in_cell;

% desired order (if any SVs are added, they should also be inserted here)
for aa = 1:num_SV_types
    if strcmp(rank_cell{aa,2},'power')==1 && strcmp(rank_cell{aa,3},'abs_diff')== 1
        rank_cell{aa,1} = 1;
    elseif strcmp(rank_cell{aa,2},'power')==1 && strcmp(rank_cell{aa,3},'abs_sum')== 1
        rank_cell{aa,1} = 2;
    elseif strcmp(rank_cell{aa,2},'log_power')==1 && strcmp(rank_cell{aa,3},'abs_sum')== 1
        rank_cell{aa,1} = 3;
    elseif strcmp(rank_cell{aa,2},'log_power')==1 && strcmp(rank_cell{aa,3},'abs_diff')== 1
        rank_cell{aa,1} = 4;    
    elseif strcmp(rank_cell{aa,2},'croot_power')==1 && strcmp(rank_cell{aa,3},'abs_sum')== 1
        rank_cell{aa,1} = 5;
    elseif strcmp(rank_cell{aa,2},'croot_power')==1 && strcmp(rank_cell{aa,3},'abs_diff')== 1
        rank_cell{aa,1} = 6;       

    elseif strcmp(rank_cell{aa,2},'azimuth_probabilities')==1 && strcmp(rank_cell{aa,3},'pearsons_corr')== 1
        rank_cell{aa,1} = 7;
    elseif strcmp(rank_cell{aa,2},'log_azimuth_probabilities')==1 && strcmp(rank_cell{aa,3},'pearsons_corr')== 1
        rank_cell{aa,1} = 8;
    elseif strcmp(rank_cell{aa,2},'azimuth_probabilities')==1 && strcmp(rank_cell{aa,3},'pearsons_corr_context')== 1
        rank_cell{aa,1} = 9;
    elseif strcmp(rank_cell{aa,2},'log_azimuth_probabilities')==1 && strcmp(rank_cell{aa,3},'pearsons_corr_context')== 1
        rank_cell{aa,1} = 10;

        
    elseif strcmp(rank_cell{aa,2},'periodicity_degree')==1 && strcmp(rank_cell{aa,3},'pearsons_corr')== 1
        rank_cell{aa,1} = 11;
    elseif strcmp(rank_cell{aa,2},'periodicity_degree')==1 && strcmp(rank_cell{aa,3},'general_corr')== 1
        rank_cell{aa,1} = 12;
    elseif strcmp(rank_cell{aa,2},'normalized_ac')==1 && strcmp(rank_cell{aa,3},'pearsons_corr')== 1
        rank_cell{aa,1} = 13;
    elseif strcmp(rank_cell{aa,2},'normalized_ac')==1 && strcmp(rank_cell{aa,3},'general_corr')== 1
        rank_cell{aa,1} = 14;
    elseif strcmp(rank_cell{aa,2},'periodicity_degree')==1 && strcmp(rank_cell{aa,3},'pearsons_corr_context')== 1
        rank_cell{aa,1} = 15;
    elseif strcmp(rank_cell{aa,2},'periodicity_degree')==1 && strcmp(rank_cell{aa,3},'general_corr_context')== 1
        rank_cell{aa,1} = 16;     
    elseif strcmp(rank_cell{aa,2},'normalized_ac')==1 && strcmp(rank_cell{aa,3},'mean_max')== 1
        rank_cell{aa,1} = 17;
    elseif strcmp(rank_cell{aa,2},'time_sig')==1 && strcmp(rank_cell{aa,3},'general_corr')== 1
        rank_cell{aa,1} = 18;
    elseif strcmp(rank_cell{aa,2},'normalized_ac')==1 && strcmp(rank_cell{aa,3},'pearsons_corr_context')== 1
        rank_cell{aa,1} = 19;
    elseif strcmp(rank_cell{aa,2},'normalized_ac')==1 && strcmp(rank_cell{aa,3},'general_corr_context')== 1
        rank_cell{aa,1} = 20;           
    
    else
        error('unknown SV type')
    end     
end 
                                
ranked_cell         = sortrows(rank_cell,1);
out_cell            = ranked_cell(:,2:3);

end

