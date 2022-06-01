function [selection,remainder] = drawfrompool(pool,num_draw)
% "ziehen ohne Zurücklegen" (?)
%Trekken zonder teruglegging

if nargin < 2
    num_draw = 1;
end

if num_draw > numel(pool)
    error('Number of draws exceeds pool size')
end

selection = zeros(1,num_draw);

for i = 1:num_draw
    n = numel(pool);
    ind = randi(n);
    draw = pool(ind);
    pool(ind) = [];
    selection(i) = draw;
end

remainder = pool;