function erb_filterbank = create_erb_filterbank(Fs, nfft, num_filters, f_low, f_high)
% create a filterbank in matrix form
% for each filter the matrix contains one column
% the filters are triangular and have a 50% overlap
%
% input nfft = length of fft output
%       Fs = sampling frequency in Hz
%       num_filters = number of filters
%       f_low = lower cutoff frequency of the filterbank
%       f_high = upper cutoff frequency of the filterbank
%
% width of bands = 1 ERB (make this an input argument?)
% 
% the ERB bands should overlap (triangular window)
% construct frequency vector (assign Hz values to the bins in the fft output)
% spacing between the frequency samples/ bins is fs/nfft (the fft size)
%
% this function uses other functions:
%       - inverbrate.m
%       - Erbw.m
%       - erbrate.m

f_min = 0;  % lowest frequency in fft output
f_max = 0.5*Fs; % highest freq in fft output
K = nfft/2 + 1; % length of each filter (half of the fft output)
f = linspace(f_min, f_max, K); % vector with frequencies in Hz that correspond to the fft output

cutoffs_erb = erbrate(f_low)+[0:(num_filters+1)]*((erbrate(f_high)-erbrate(f_low))/(num_filters+1)); % vector with cutoff freqs of filters in ERB
cutoffs_hz = inverbrate(cutoffs_erb);

% make sure that each filter has at least 2 ERB width (1 ERB overlap)
erbstep = cutoffs_erb(4)-cutoffs_erb(3);
if erbstep <= 1 % erbstep should be at least 1 to give one filter a bandwidth of 2ERB
    error('Bandwidth of erb filters is smaller than 2. Decrease number of filters.')
end

erb_filterbank = zeros(num_filters,K); % length of each filter x number of filters

for m = 1:num_filters
    
    for k = 1:K; 
        
        f_k = f(k);
        
        if (f_k >= (cutoffs_hz(m))) && (f_k <= cutoffs_hz(m+1))
            erb_filterbank(m,k) = (f_k-cutoffs_hz(m))/(cutoffs_hz(m+1)-cutoffs_hz(m)); 
            
        elseif (f_k >= cutoffs_hz(m+1)) && (f_k <= cutoffs_hz(m+2))
            erb_filterbank(m,k) = (cutoffs_hz(m+2)-f_k)/(cutoffs_hz(m+2)-cutoffs_hz(m+1));        
        else
            erb_filterbank(m,k) = 0;
        end
    end
    
    erb_filterbank(m,:)=erb_filterbank(m,:).*(1/sum(erb_filterbank(m,:))); % normalize the filters
end



end




