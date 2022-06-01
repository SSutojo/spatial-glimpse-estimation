function [gtfb, filtStatOut] = OZGTFB(fs, fCenter, signal, filtStatIn)
%OZGTFB.m One Zero Gammatone Filter Bank implemenation.
%
%> This function implements a One-Zero-All-Pole Gammaton filter bank based
%> on the implementation of Pflueger and Hoeldrich, 1997.
%> Each gammatone filter consists of four second order IIR-filters except
%> the very first stage. As shown in Pflueger1997 the addition one zero
%> there has a high improves of the filter performance near DC.
%
%> @param fs Sampling frequency.
%> @param fCenter Center frequencies of the gammatone filters.
%> @param signal Audio input signal which has only one channel.
%> @param filtStatIn Filter states for block-wise calculation.
%> @param bAlign Enables/Disables peak alignment at the end of calculation
%> @retval gtfb Gammatone filter outputs for all channels. 
%> @retval filtStatOut Filter states at the end of calculation.

% Check input argument bAlign.
if nargin < 5
    bAlign = 0;
end

% Static paramter.
T         = 1/fs;
earQ      = 9.26449;           
bwMin     = 24.7;
bwCorrect = 1.019;
nSamp     = size(signal, 1);

nFilt = length(fCenter);

% Preallocate variables.
filtOut     = zeros(nSamp, nFilt);
filtStatOut = zeros(size(filtStatIn));
% Calculate GTFB paramter.
fCentErb = ((fCenter/earQ) + bwMin);
B        = 2 * pi * fCentErb * bwCorrect;
twoPiT   = 2 * pi * T;
z        = exp(1i * twoPiT * fCenter);

% Calculate center frequency gain for two-pole filters and first stage gain with added zero. 
threeStageGain = abs((1 - 2 * cos(twoPiT * fCenter) .* z .* exp(-B * T) + exp(-2 * B *T) .* z.^2));    
firstGain      = abs(threeStageGain ./ (1 - z));
gain           = firstGain .* threeStageGain.^3;

% Generate forward and feedback filter coefficent.
% SOS:    H(Z) = (b0 + b1 * z) / (a0 + a1 * z + a2 * z^2);
% Coeff:  b0 = 1; b1 = {-1, 0}; a0 = 1; a1 = -2cos(2 pi T fc)e^(-BT); a2 = e^(-2BT); 
% Note that b1 = -1 only for first filter stage and zero otherwise.
forward       = 1;
oneZero       = -1;
feedback      = zeros(nFilt, 3);
feedback(:,1) = 1; 
feedback(:,2) = -2 * cos(twoPiT * fCenter) .* exp(-B * T);
feedback(:,3) = exp(-2 * B * T);

% Filter signal with four second-order stages.
for iFilt = 1:nFilt
    [stage1, filtStatOut(:, 1, iFilt)] = filter([forward, oneZero], feedback(iFilt, :), signal, filtStatIn(:, 1, iFilt));
    [stage2, filtStatOut(:, 2, iFilt)] = filter(forward, feedback(iFilt, :), stage1, filtStatIn(:, 2, iFilt));
    [stage3, filtStatOut(:, 3, iFilt)] = filter(forward, feedback(iFilt, :), stage2, filtStatIn(:, 3, iFilt));
    [filtOut(:, iFilt), filtStatOut(:, 4, iFilt)] = filter(forward, feedback(iFilt, :), stage3, filtStatIn(:, 4, iFilt));
end

% Normalize filter gain to 0dB at center frequency.
gtfb = filtOut .* gain(ones(nSamp, 1), :);
