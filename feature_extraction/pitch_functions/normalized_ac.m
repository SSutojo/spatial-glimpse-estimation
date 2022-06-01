function NAC = normalized_ac(sig_in, period_start, period_end, period_step,Fs)
% calculate the normalized autocorrelation function for the signal sig_in.
% The time lag is varied within the interval from lag_start to lag_end

% INPUTS : sig_in: time signal
%          lag_start, lag_end : interval in which time lag is varied
% OUTPUTS: NAC: a vector that contains the normalized autocorrelation for
% each of te L values (time lag values)

n_period_start = round(period_start*10^(-3)*Fs); % change the ms values for period start/end/step into sample numbers
n_period_end   = round(period_end  *10^(-3)*Fs); 
n_period_step  = round(period_step *10^(-3)*Fs);

upper_bound = length(sig_in)-1; % change this
period_vec = n_period_start:n_period_step:n_period_end;
NAC = zeros(length(period_vec),1);

for idx = 1:length(period_vec); % loop over the period candidates
    L = period_vec(idx);
    n_prod = [];
    sig_sq = [];
    delayed_sq = [];
    
    for j = 1:(upper_bound-L);              % loop over all samples, excluding L samples at the end of sig_in
        n_prod(j) = sig_in(j)*sig_in(j+L);  % product of the original signal and its delayed version
        sig_sq(j) = sig_in(j).^2;           % squared original signal
        delayed_sq(j) = sig_in(j+L).^2;     % squared delayed signal
    end
    
    NAC(idx) = sum(n_prod)/(sqrt(sum(sig_sq)+eps)*sqrt(sum(delayed_sq))); % normalized AC for each L value  PROBLEM sig_sq can be 0, to avoid this add eps
    
end

end