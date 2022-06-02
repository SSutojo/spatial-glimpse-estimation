function center_freqs = ERB_center_freqs(f_low, f_high, erbstep)
% create a vector with center frequencies of ERB filters between the lowest
% given frequency in Hz (f_low) and the highest given frequency (f_high).
% The ERB filters are spaced according to specified erbstep.
%
% inputs: f_low : lowest frequency in Hz
%         f_high: highest frequency in Hz
%         erbstep: number of ERBs in one filter
% outputs: center_freqs : a vector with center frequencies in Hz


HzToERBRate         = @(Hz) 21.4*log10(Hz.*0.00437 + 1); 
ERBRateToHz         = @(erbf) (10.^(erbf/21.4)-1)/4.37e-3;

er_low              = ceil(HzToERBRate(f_low));  % determine the erbrate of the lowest given frequency
er_low              = er_low + (erbstep/2);      % increase er_low by half a band, so that the lowest filter does not extend below f_low

er_high             = floor(HzToERBRate(f_high))- (erbstep/2); % culculate the highest erb and decrease by half an erbstep so that the filter won't exceed the highest frequency

erb_center_freqs    = [er_low:erbstep:er_high];
center_freqs        = ERBRateToHz(erb_center_freqs); 


