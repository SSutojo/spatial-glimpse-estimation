function [PD, NAC, CFR] = periodicity_feature(sig_in, period_start, period_end, period_step, Fs)
% calculate the periodicity feature to combines the normalized
% autocorrelation (NAC) and the comb filter ratio (CFR)
% according to:
% Z. Chen and V. Hohmann, "Online Monaural Speech Enhancement Based on
% Periodicity Analysis and A Priori SNR Estimation", IEEE Trans. Audio. 
% Speech, Lang. Process., vol. 23, no. 11, pp.1904-1916, 2015.
%
% 
% INPUTS: period_start: smallest period/ time lag (in ms)
%         period_end : longest period in ms
%         period_step: intervals between the tested period candidates
%         sig_in : time signal
% OUTPUTS: 
%         PD:   periodicity degree
%         NAC: normalized autocorrelation
%         CFR: comb-filter-ratio
%

NAC = normalized_ac(sig_in, period_start, period_end, period_step,Fs); % a vector containing all the NAC values for period candidates between lag_start and lag_end
CFR = comb_filter_ratio(sig_in, period_start, period_end, period_step,Fs); % all the CFR values


PD = max(0.01 ,NAC.*CFR); % calculation of the actual periodicity feature


end