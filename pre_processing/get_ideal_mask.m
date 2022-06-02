function [source_map, IBMs, local_SNRs] = get_ideal_mask(power_mat_sep, opt_mask, bin_mask)
%calculate the ideal source map and ideal binary masks for each source by
%comparing the single sources' powers in each time-frequency unit  
%
% INPUTS: - power_mat_sep : (number of time frames X number of Filterbands X number of sources)
%           energy per t-f-unit for each component
%           power_mat_sep has only one ear Channel, summed over both ears or monaural
%         - opt_mask: option ('y' or 'n') if a binary mask (e.g. a VAD) should be applied or not
%         - bin_mask : the binary mask that is optionally applied. size(num. time frames X num. filterbands)
%         
% OUTPUTS:- source_map: matrix that indicates which source (signal component)
%           is dominant in each time-frequency unit (size: time frames x filterbands)
%           the ID given each speaker is the same as its row number in the input
%           matrix "signal_components"
%         - IBMs: ideal binary masks for each of the signal components. The
%           IBM is set to 1 if the corresponding source is dominating whithin
%           the specified t-f-unit. Size (num. time frames X num. filterbands X num. sources)
%         - local_SNRs: has the same dimensions as IBMs matrix but instead
%           of 0s and 1s contains the SNR for the regarded source in each
%           t-f-unit


num_sources = size(power_mat_sep,3);                % get the number of signal components in the mixture

% check if the VAD mask should be applied
if opt_mask == 'y'                                    
    max_val = max(max(max(power_mat_sep)));         % get the maximum value in power_mat
    VAD_mask2 = abs((bin_mask-1)*(max_val+1));      % the active t-f-units are set to zero, the silent parts are set to (max_val+1)  
    power_mat_sep_masked = cat(3,power_mat_sep,VAD_mask2); % add the mask as another source. 
                                                    % The source map will be set to the ID of the mask in units where the bin_mask is zero    
elseif opt_mask == 'n'
    power_mat_sep_masked = power_mat_sep;                  % leave the power matrix as it is
end


% compare the sources' powers in each t-f-unit to get the dominant source and its local SNRs
[~, source_map] = max(power_mat_sep_masked(:,:,:),[],3);

IBMs = zeros(size(source_map,1),size(source_map,2),num_sources);
local_SNRs = zeros(size(IBMs));
for ii = 1:num_sources
    IBMs(:,:,ii)= source_map == ii;
    noise_mat   = power_mat_sep(:,:,setdiff(1:end,ii));
    noise_sum   = sum(noise_mat,3);
    local_SNRs(:,:,ii) = 10*log10(power_mat_sep(:,:,ii)./noise_sum);
end


end