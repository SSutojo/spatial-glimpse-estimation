function [ind, val] = findNearest(x, y)
%findNearest   Returns indices of values y that are nearest to the values 
%   in vector x.
%
%USAGE
%   [IDX,VAL] = findNearest(X,Y);
%
%INPUT ARGUMENTS
%     X : list of values
%     Y : values which should be found in X
%
%OUTPUT ARGUMENTS
%   IDX : indices of nearest values of Y in X.
%   VAL : nearest values of Y in X.
%
%EXAMPLE
%   % Create value list
%   ref = 0:1e-3:1;
% 
%   % Find closest values of 3 randomized numbers in ref
%   findNearest(ref,rand(3,1))

x=x(:);     % make column
y=y(:).';   % make row

[val, ind] = min(abs(y(ones(size(x,1),1),:)-x(:,ones(1,size(y,2)))));

% Get values
if nargout > 1
    val = x(ind);
end