function out = mono(sig)
%mono   Mono downmix of the audio object.

%   Author  :   Tobias May (tobias.may@uni-oldenburg.de)
%
%   Version :   0.1     2008/05/07 
%   ***********************************************************************

% Check if IN is an audio structure
if ~isAudio(sig);
   error('Input argument IN must be an audio structure.') 
end
    
% Copy signal structure
out           = sig;
% Mono downmix
out.data      = mean(sig.data,2);
out.nChannels = 1;


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************