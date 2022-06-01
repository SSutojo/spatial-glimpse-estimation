function [bacc] = evaluate_IBM(true_IBM,est_IBM)
%EVALUATE_IBM calculate the overlap/ balanced accuracy between the true
%ideal binary mask (IBM) and the estimated IBM
% INPUTS
%       true_IBM    : ground truth IBM of a source. Dim(num_win X num_bands)
%       est_IBM     : estimated IBM (same dimensions as true_IBM)
% OUTPUTS
%       bacc        : balanced accuracy

if size(true_IBM)~=size(est_IBM)
    error('true IBM and est_IBM must have the same dimensions')
end

P = sum(sum(true_IBM)); % total number of Positioves = 1s (t-f-units labeled with 1) in the true IBM
N = numel(true_IBM)-P; % total number of Negatives = zeros in true IBM

overlap_mat = true_IBM.*2-est_IBM;

TP = length(find(overlap_mat==1));  % true positives  = correctly detected ones
FP = length(find(overlap_mat==-1)); % false positives = falsely detected ones
FN = length(find(overlap_mat==2));  % false negatives = undetected ones
TN = numel(true_IBM)-(TP+FP+FN);    % true negatives  = correctly detected zeros

bacc = (TP/P + TN/N)/2;

end

