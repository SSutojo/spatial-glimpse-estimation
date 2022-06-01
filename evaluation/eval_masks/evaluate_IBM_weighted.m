function [weighted_bacc] = evaluate_IBM_weighted(true_IBM,est_IBM, target_power_mat, masker_power_mat, SNR_criterion)
%EVALUATE_IBM_WEIGHTED calculate the overlap/ balanced accuracy between the true
%ideal binary mask (IBM) and the estimated IBM. The t-f-units are weighted
%according to their local SNR. Some
%t-f-units contain very few overall power and are therefore (more or less
%randomly) labeled with 1 because there is little masker energy as
%well. Since these t-f-units contain little reliable features or target energy they 
%are less important/weighted for the IBM
% INPUTS
%       true_IBM    : ground truth IBM of a source. Dim(num_win X num_bands)
%       est_IBM     : estimated IBM (same dimensions as true_IBM)
%       target_power_mat : signal power per t-f-unit. Dim(num_win X num_bands)
%       masker_power_mat : signal power per t-f-unit of masker. Same dimensions 
%       SNR_criterion : criterion above which the t-f-units are assigned a higher weight
% OUTPUTS
%       weighted_bacc        : balanced accuracy

if size(true_IBM)~=size(est_IBM)
    error('true IBM and est_IBM must have the same dimensions')
end

P = sum(sum(true_IBM)); % total number of Positioves = 1s (t-f-units labeled with 1) in the true IBM
N = numel(true_IBM)-P; % total number of Negatives = zeros in true IBM

overlap_mat = true_IBM.*2-est_IBM;

TP = length(find(overlap_mat==1));  % true positives  = correctly detected ones
FP = length(find(overlap_mat==-1)); % false positives = falsely detected ones
FN = length(find(overlap_mat==2));  % false negatives = undetected ones
TN = length(find(overlap_mat==0));  % true negatives  = correctly detected zeros

local_SNRs = target_power_mat./masker_power_mat;

low_SNR_TP = length(find(overlap_mat==1 & local_SNRs<SNR_criterion));  % true positives with a low local SNR (get a low weight cause less important)
high_SNR_TP= length(find(overlap_mat==1 & local_SNRs>=SNR_criterion)); % true postives that get a higher weight cause they contain more information
low_SNR_TN = length(find(overlap_mat==0 & local_SNRs<SNR_criterion));
high_SNR_TN= length(find(overlap_mat==0 & local_SNRs>=SNR_criterion));

low_SNR_P  = length(find(true_IBM==1 & local_SNRs<SNR_criterion));      % positives (ones in true IBM) with low SNR
high_SNR_P = length(find(true_IBM==1 & local_SNRs>=SNR_criterion));
low_SNR_N  = length(find(true_IBM==0 & local_SNRs<SNR_criterion));
high_SNR_N = length(find(true_IBM==0 & local_SNRs>=SNR_criterion));

if low_SNR_P==0 
    low_SNR_P = 1;
elseif high_SNR_P==0
    high_SNR_P = 1;
elseif low_SNR_N==0
    low_SNR_N = 1;
elseif high_SNR_N==0
    high_SNR_N = 1;
end



weighted_bacc = ((low_SNR_TP/low_SNR_P) + (high_SNR_TP/high_SNR_P)*2 +(low_SNR_TN/low_SNR_N) + (high_SNR_TN/high_SNR_N)*2)/6;

% check if the calculation is right
if (low_SNR_TP+high_SNR_TP)~=TP
    error(':-(')
end
if (low_SNR_TN+high_SNR_TN)~=TN
    error(':-(')
end

end

