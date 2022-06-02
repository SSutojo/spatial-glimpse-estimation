function bw=Erbw(f)
% calculate equivalent rectangular bandwidth according to Glasberg and Moore(1990)
% function bw=Erbw(f);
%
% IN    f: center frequency of the filter in Hz		
% OUT   bw: bandwidth of filter in Hz

error(nargchk(1,1,nargin));

bw=(f.*0.00437+1).*24.7;
