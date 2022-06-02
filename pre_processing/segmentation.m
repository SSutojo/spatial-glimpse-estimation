function [sig_seg, num_win, cut_sig_in] = segmentation(sig_in, Fs, Lwin, hopsize)
% input signal sig_in is segmented into time frames
% 
% inputs: sig_in - time signal
%         Fs - sampling frequency
%         Lwin - desired length of one time frame in ms
%         noverlap - desired overlap of the time frames (usually here 0.5 Lwin)
%
% output: sig_seg - matrix with segmented signal (each column hold a windowed segment of sig_in)

sig_in          = sig_in(:); % sig_in should be a column vector

Lwin            = Lwin*10^(-3)*Fs;         % desired length of time frames in samples
nhopsize        = fix(hopsize*Lwin);   % number of overlapping samples
num_win         = floor((length(sig_in)-(Lwin-nhopsize))/nhopsize); % number of time frames, discard the last few samples
final_length    = num_win*nhopsize + (Lwin-nhopsize);

%padded_sig_in   = [sig_in ; zeros((final_length-length(sig_in)),1)];  % zero pad the input signal so that it can be evenly divided by Lwin
cut_sig_in   = sig_in(1:final_length,:);  % cut the input signal so that it can be evenly divided by Lwin

sig_seg = zeros(Lwin,num_win);  % mat to hold segments of sig_in
for idx = 1:num_win
    segment         = cut_sig_in(fix(1+(idx-1)*nhopsize) : (fix((idx-1)*nhopsize)+Lwin));
    sig_seg(:,idx)  = segment;
end


end