function er=erbrate(f)
% calculate erbrate according to Glasberg and Moore (1990)
%
% IN              f :   center frequency of the filter.  		
% OUT            er :   erbrate
%            

error(nargchk(1,1,nargin));

er=log(f.*0.00437+1)./0.1079;
