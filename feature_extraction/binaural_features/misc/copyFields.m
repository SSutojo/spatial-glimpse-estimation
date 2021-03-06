function out = copyFields(out,in)

%   Developed with Matlab 7.3.0.267 (R2006b). Please send bug reports to:
%   
%   Author  :   Tobias May  tobias.may@uni-oldenburg.de
%
%   Version :   1.0   2007/08/17
% *************************************************************************

if nargin < 2 || isempty(in); help(mfilename); return; end

% Get all fieldnames of struct 'in'
inFields = fieldnames(in);

% Loop over all fieldnames 
for ii = 1 : length(inFields)
    % Copy field if the fieldname also exists in the structure 'out'
    if isfield(out,inFields{ii})
        out.(inFields{ii}) = in.(inFields{ii});
    end
end