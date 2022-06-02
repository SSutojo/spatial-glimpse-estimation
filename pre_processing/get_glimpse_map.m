function [glimpse_map] = get_glimpse_map(regions_est)
%GET_GLIMPSE_MAP Summary of this function goes here
%   Detailed explanation goes here

num_win = ceil(size(regions_est,1)/2);
nFilts  = ceil(size(regions_est,2)/2);

glimpse_map             = zeros(num_win,nFilts);
for bb = 1:num_win % time frames
    for cc = 1:nFilts % filterbands
        glimpse_map(bb,cc) = regions_est((2*bb)-1,(2*cc)-1);
    end
end
end

