function CFR = comb_filter_ratio(sig_in, period_start, period_end, period_step,Fs)
% calculate the comb filter ratio (CFR)
% CFR = ratio of the frame energy of the summation of the original signal and its delayed version to the frame energy
% of the difference between original signal and the delayed version
% the CFR has a peak at the delay corresponding to P_o (the period
% corresponding to the fundamental frequency)
% INPUTS : sig_in : signal that shall be analyzed
%          lag_start and lag_end : mark the interval in which the time lag
%          should be varied to find the value P_o at which the CFR exhibits
%          a peak
% OUTPUTS: cfr : comb filter ratio

n_period_start = round(period_start*10^(-3)*Fs); % change the ms values for period start/end/step into sample numbers
n_period_end   = round(period_end  *10^(-3)*Fs); 
n_period_step  = round(period_step *10^(-3)*Fs);

upper_bound = length(sig_in)-1;
period_vec = n_period_start:n_period_step:n_period_end;
CFR = zeros(length(period_vec),1);
for idx = 1:length(period_vec);
    L = period_vec(idx);
    n_sum = [];
    n_dif = [];
    for j = 1:(upper_bound-L);
        n_sum(j) = sig_in(j)+sig_in(j+L); % sum of original signal and its delayed version
        n_dif(j) = sig_in(j)-sig_in(j+L); % difference between original signal and delayed version
    end
    CFR(idx) = sum(n_sum.^2)/sum(n_dif.^2);
    
end 



end


