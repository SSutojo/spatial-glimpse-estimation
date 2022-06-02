function f=inverbrate(er)
% calculate center frequency of filter from erbrate according to Glasberg and Moore (1990)
% IN   er : erbrate corresponding to the center frequency of the filter.  		
% OUT  f : center frequency in Hz  
%

error(nargchk(1,1,nargin));

f=(exp(er.*0.1079)-1)./0.00437;
