function data = compress(data,method)
%compress   Apply compression to input signal.
%   
%USAGE
%     data = compress(data,method)
%
%INPUT ARGUMENTS
%     data : N-dimensional data
%   method : string specifying type of compression
%            'log'        - natural logarithm
%            'log_md'     - natural logarithm with an offset of 1, such
%                           that 0s do not produce INFs.
%            'log10'      - base 10 logarithm
%            'loga'       - scale-able logarithmic compression 
%            'cuberoot'   - cuberoot compression ^(1/3) 
%            'mag'        - 20log10-compression with 80 dB dynamic range
%            'pow'        - 10log10-compression with 80 dB dynamic range 
%            'nonlinear'  - non-linear hair cell compression [1]
%            '1over15'    - power law compression with ^(1/15) [1]
%            'inst'       - instantaneous compression with ^0.4 [2]
%            'htan'       - hyperbolic tangent compression [3]
%            'htan_inv'   - inverse hyperbolic tangent compression [3]
% 
%OUTPUT ARGUMENTS
%     data : compressed N-dimensional data
% 
%REFERENCES 
%   [1] N. Moritz, J. Anemüller and B. Kollmeier. "Amplitude Modulation
%       Filters as Feature Sets for Robust ASR: Constant Absolute or
%       Relative Bandwidth?", Interspeech, 2012.    
% 
%   [2] S. D. Ewert and T. Dau. "Characterizing frequency selectivity for
%       envelope fluctuations", JASA, vol. 108(3), pp. 1181-1196, 2000.
% 
%   [3] D. S. Williamson and D. Wang. "Time-Frequency Masking in the
%       Complex Domain for Speech Dereverberation and Denoising", IEEE/ACM
%       Transactions on Audio, Speech, and Language Processing, vol. 25(7),
%       pp. 1492-1501, 2017.  

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013-2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2013/08/11
%   v.0.2   2014/04/02 added non-linear compression
%   v.0.3   2017/07/06 added hyperbolic tangent compression
%   ***********************************************************************
   
    
%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% APPLY COMPRESSION
%
%
% Select method
switch(lower(method))
    case {'' 'nothing'}
        % data = data;
         
    case 'log'
        data = log(data);
        
    case 'log_md'
        % Add one, such that a lower bound of zero can be used in the
        % context of missing data classification
        data = log(data + 1);
        
    case 'log10'
        data = log10(data);
        
    case 'inst'
        % Instantaneous compression with an exponent of 0.4
        data = sign(data) .* abs(data).^0.4;
        
    case 'loga'
        a = 3;
        data = log10(1 + a * data) ./ log10(1 + a);
        
    case 'cuberoot'
        data = power(data,1/3);
    
    case 'mag'
        % Limit dynamic range to 80 dB
        dynamicRange = 80;
        
        data = 20 * log10(data + 10.^(-dynamicRange/20));
            
    case 'pow'
        % Limit dynamic range to 80 dB
        dynamicRange = 80;
        
        data = 10 * log10(data + 10.^(-dynamicRange/10));
                    
    case 'nonlinear'
        % Non-linear operation to mimic outer hair cell compression [1]
        data = (power(data,0.4) + log(data) + 1)/2;
		
    case '1over15'
        % Non-linear operation to mimic outer hair cell compression [1]
        data = power(data,1/15);		
        
    case 'htan'
        
        % Hyperbolic tangent compression
        Q = 1; C = 0.5;
        
        data = Q * (1 - exp(-C .* data)) ./ (1 + exp(-C .* data));
        
    case 'htan_inv'
        
        % Hyperbolic tangent compression
        Q = 1; C = 0.5;
        
        data = -1/C * log((Q - data) ./ (Q + data));
        
    otherwise
        error('Type of compression is not supported!')
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