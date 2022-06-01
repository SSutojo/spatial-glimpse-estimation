function F = genFilterABCD(fs,weight)
%genFilterABCD   Design weighting filter.
%
%USAGE 
%   F = genFilterABCD(fs)
%   F = genFilterABCD(fs,weight)
%
%INPUT ARGUMENTS
%       fs : sampling frequency (Hz)
%   weight : specify frequency weighting: 'A' 'B' 'C' or 'D' 
%            (default, weight = 'A')
%
%OUTPUT ARGUMENTS
%   F : filter object
%
%NOTE 
%   This function is based on weightings.m, a function of the Psysound3
%   Toolbox (see PsySound.org). 
%
%EXAMPLE
%   genFilterABCD(48e3,'A')
% 
%   See also genFilterRemoveDC.

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2009 
%              TUe Eindhoven 
%              t.may@tue.nl    
%
%   History :
%   v.0.1   2009/09/26
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
narginchk(1, 2);

% Set default values
if nargin < 2 || isempty(weight); weight = 'A'; end


%% ***************************  FILTER DESIGN   ***************************
% 
% 
% NFFT resolution for plotting 
nfft = 1024;

% The frequency of each Fourier bin
f = (0:fs/2/(nfft-1):fs/2);

% Set up the s-plane variable
s = 1i*2*pi*f; 
    
% Determine the weighting filter frequency responses.
switch upper(weight)
    case 'A' % A-weighting filter
        K = 7.39705e9;
        freqResp = K*s.^4./((s+129.4).^2 .* (s+676.7).*(s+4636)...
                   .*(s+76655).^2);
        
        zrs =  [0; 0; 0; 0];
        pls = -[129.4; 129.4; 676.7; 4636; 76655; 76655];
        
    case 'B' % B-weighting filter
        K = 5.99185e9;
        freqResp = K*s.^3./((s+129.4).^2 .* (s+995.9).*(s+76655).^2);
        
        zrs =  [0; 0; 0];
        pls = -[129.4; 129.4; 995.9; 76655; 76655];
        
    case 'C' % C-weighting filter
        K = 5.91797e9;
        freqResp = K*s.^2./((s+129.4).^2 .*(s+76655).^2);
                
        zrs =  [0; 0];
        pls = -[129.4; 129.4; 76655; 76655];
        
    case 'D' % D-weighting filter
        K = 91104.32;
        freqResp = K*s.*(s.^2+6532*s+4.0975e7)./((s+1776.3) .* ...
                   (s+7288.5).*(s.^2+21514*s+3.8836e8));
        
        zrs = [0; roots([1 6532 4.0975e7])];
        pls = [-1776.3; -7288.5; roots([1 21514 3.8836e8])];
        
    otherwise % Unknown request
        error(['Filter weight ''',weight,''' is unknown. ',...
               'Options are ''A'', ''B'', ''C'' or ''D'''])
end

% Generate the filter
if (2*f ~= fs) 
    % Correct small frequency error on the last fourier sample.
    f(end) = fs/2;
end

% Use the bilinear transformation to discretize the transfer function. 
warning('off','MATLAB:nearlySingularMatrix');
[Zd, Pd, Kd] = bilinear(zrs, pls, K, fs);
warning('on','MATLAB:nearlySingularMatrix');
[b, a]       = zp2tf(Zd, Pd, Kd);


%% **********************  CREATE FILTER STRUCTURE   **********************
% 
% 
% Create filter object
F = genFilterObj(b,a,'fs',fs,'label',[upper(weight),'-weighting filter']);


%% ************************  SHOW FILTER RESPONSE   ***********************
% 
% 
% Show result
if nargout == 0
    % Plot frequency respones
    mag = freqz(b, a, 2*pi*f/fs);
    
    figure;
    h = semilogx(f, 10*log10(eps + abs([mag; freqResp]).^2));
    set(h(1),'LineStyle','-','Color','k','LineWidth',2)
    set(h(2),'LineStyle','--','Color',[0 0.5 0],'LineWidth',2)
    grid on;
    xlim([0 fs/2])
    ylim([-50 20])
    title([weight, ' weighted']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend('original','filter response', 'Location', 'NorthWest');
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