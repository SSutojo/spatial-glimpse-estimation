function [locMaxIdx,maxIdx,locMaxIter,maxIter] = findLocalPeaks(in,mode)
%findLocalPeaks   MEX-File for local and global maximum search.
% 
%   A value is considered as local/global peak if the adjacent amplitudes 
%   decrease.
% 
%USAGE
%   [LOCMAXIDX,MAXIDX,LOCMAXITER,MAXITER] = findLocalPeaks(IN,MODE);
% 
%INPUT ARGUMENTS
%     IN : input matrix [nSamples x nChannels]
% 	MODE : mode for maxima search
% 	       1 = endpoints are not accepted as local or global maxima. So the
% 	       number of global maxima does not necessarily corresponds to the 
%          number of channels as there might be channels without any 
%          maximum.
%
% 	       2 = also endpoints are accepted as local or global maxima. The 
%          number of global maxima always corresponds to the number of 
%          channels. (default, MODE = 2)
% 
%OUTPUT ARGUMENTS
%   LOCMAXIDX : indices of all local maxima. If IN is a matrix LOCMAXIDX 
%               can be used to directly access the maxima.
% 	   MAXIDX : index of global maximum per channel
%  LOCMAXITER : indices identifying the channel of local maxima.
%     MAXITER : indices identifying the channel of global maxima.
% 
%EXAMPLE
%   % Create input data
%   blockSize = 100;
% 	in        = rand(blockSize,10);    
% 	in        = filter([0.15 0.15],[1 -0.6],in);
% 
% 	% Find local and global maxima
% 	[a,b,c,d] = findLocalPeaks(in);
% 	% Use local and global indices for plotting on a 2D surface
% 	a = mod(a-1,blockSize)+1;
% 	b = mod(b-1,blockSize)+1;
% 
% 	% Plot 2D surface
% 	figure;hold on;
% 	imagesc(in);axis xy
% 	plot(c,a,'xw','MarkerSize',4)
% 	plot(d,b,'*w','MarkerSize',8)
% 	axis tight; colorbar; 

%   Developed with Matlab 7.5.0.342 (R2007b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2008 
%              TUe Eindhoven 
%              t.may@tue.nl    
%
%   History :
%   v.0.1   2008/07/14
%   ***********************************************************************

% Check for proper input arguments
narginchk(1,2);

% Set default values...
if nargin < 2 || isempty(mode);  mode = 2; end

% TODO % Handle plateaus, which are currently not considered as peaks

% Find local and global peaks (MEX processing)
[locMaxIdx,maxIdx,locMaxIter,maxIter] = findLocalPeaksMEX(in,mode);


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