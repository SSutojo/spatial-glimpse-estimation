function c=logtimes(a,b)
% function c=logtimes(loga,b) returns log(a*b) and a is passed
% as log(a). useful for probabilistic calculations when a is a density
% and can be very small.

%       Copyright (c) Alexander G Dimitrov 2001


Maxa=max(a,[],2); a=a-repmat(Maxa,1,size(a,2));  % rescale (shift
                                                 % logs), per
                                                 % observation 
c=log(exp(a)*b)+repmat(Maxa,1,size(b,2)); % exp(a).*exp(Maxa) gets
                                          % to the old exp(a)