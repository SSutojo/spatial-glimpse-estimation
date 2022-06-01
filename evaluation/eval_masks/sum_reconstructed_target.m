function [rec_target] = sum_reconstructed_target(ideal_source_map, estimated_source_map, mixed_power_mat)
%SUM_MISMATCHED_ENERGY sum the target energy that is reconstructed when
%applying the estimated binary mask to the mixed signal, energy that's
%correctly assigned is divided by total energy
%   INPUTS:
%   ideal_source_map (#frames X #subbands)
%   estimated_source_map (same size/dims as IBM)
%   mixed_power_mat - power per T-F unit summed over both ear channels
%                     (also #frames X # subbands)
%
%   OUTPUTS:
%   rec_target - reconstructed target energy over total energy


%% get all correctly assigned T-F units and sum their power
power_mat_rec = mixed_power_mat;

power_mat_rec(find(ideal_source_map ~= estimated_source_map)) = 0;


rec_target = 10*log10(sum(power_mat_rec(:))/sum(mixed_power_mat(:))); % reconstructed target energy over total energy


end

