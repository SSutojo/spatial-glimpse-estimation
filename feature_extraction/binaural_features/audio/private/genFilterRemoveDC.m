function F = genFilterRemoveDC(fs,cutoffHz)
%genFilterRemoveDC   Design 4th order high-pass to remove DC components.
%
%USAGE 
%   F = genFilterRemoveDC(fs)
%   F = genFilterRemoveDC(fs,cutoffHz)
%
%INPUT ARGUMENTS
%         fs : sampling frequency (Hz)
%   cutoffHz : cutoff frequency (default, cutoffHz = 20)
%
%OUTPUT ARGUMENTS
%   F : filter object
%
%REFERENCES
%   [1] ITU-R Recommendation BS.1387: Method for objective measurements of
%       perceived audio quality (PEAQ)
% 
%   [2] Kabal, P. (2002) "An Examination and Interpretation of ITU-R
%       BS.1387: Perceptual Evaluation of Audio Quality". TSP Lab Technical
%       Report, updated 2003.  
%
%EXAMPLE
%   genFilterRemoveDC(44.1e3);
% 
%   See also genFilterABCD.

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2009 
%              TUe Eindhoven
%              t.may@tue.nl   
%
%   History :
%   v.0.1   2009/09/26
%   v.0.2   2009/10/07 changed filter design according to [1]
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
narginchk(1, 2);

% Set default value
if nargin < 2 || isempty(cutoffHz); cutoffHz = 20; end


%% ***********************  CREATE HIGHPASS FILTER  ***********************
% 
% 
% Note: There is a more efficient implementation based on two 2nd order
% filters using second order sections [2].

% Design 4th order butterworth highpass
[b, a] = butter(4, cutoffHz / (0.5 * fs),'high');

% Create filter object
F = genFilterObj(b,a,'fs',fs,'label','DC removal filter');


%% ***********************  SHOW FILTER RESPONSE  *************************
% 
% 
% Show filter characteristic
if nargout == 0
    plotFilter(F);
end


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