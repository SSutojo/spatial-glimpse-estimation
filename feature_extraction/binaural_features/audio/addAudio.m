function out = addAudio(varargin)
%addAudio   Add audio objects.
%   
%USAGE
%    aObj = addAudio(varargin)
%    aObj = addAudio(IN1,IN2,...)
%
%INPUT ARGUMENTS
%    varargin : audio signals arranged as:
%               1. multiple input arguments IN1,IN2,...  or 
%               2. as higher dimensional audio object (see EXAMPLE)
%
%OUTPUT ARGUMENTS
%    aObj : audio structure
%
%NOTE
%   All audio objects are simply added and the amplitude of the output
%   signal is NOT scaled by the number of signals.
%  
%   See also isAudio, updateAudio, readAudio, normalizeAudio, plotAudio 
%            and resampleAudio.
% 
%EXAMPLE
%   % Load audio material
%   load chirp
% 
%   % Create higher dimensional audio object
%   in(1) = genAudio(y,Fs,'chirp');
%   in(2) = genAudio(y,Fs,'chirp');
%   in(3) = genAudio(y,Fs,'chirp');
% 
%   % Add all three signals
%   addAudio(in);


%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
% 
%   Author  :  Tobias May, � 2007-2010 
%              University of Oldenburg and TU/e Eindhoven 
%              tobias.may@uni-oldenburg.de   t.may@tue.nl
%
%   History : 
%   v.0.1   2007/08/17
%   v.0.2   2008/05/07 allow arbitrary number of input signals
%   v.0.3   2008/11/04 support different data range between audio signals
%   v.0.4   2009/05/15 added single cell support
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
narginchk(1, inf);

try
    % Extract audio structures 
    if all(size(varargin) == 1) && isstruct(varargin{1})
        in = varargin{1};
        varargin = cell(length(in),1);
        % Convert cells to structs
        for ii = 1 : length(in)
            varargin{ii} = in(ii);
        end
    elseif size(varargin,1) == 1 && size(varargin,2) > 1 ...
                                 && isstruct(varargin{1})
        % Do nothing ...
    else
        % Support single cell containing multiple audio structures
        varargin = varargin{:};
    end
catch ME
    % Report error
    error('Supply audio structure either in a struct arrays or in a cell.')
end

% Get number of input signals
nSignals = length(varargin);
% Keep track of sampling rates
fsIdx = zeros(nSignals,1);
% Keep track of the number of channels
nChan = zeros(nSignals,1);
% Init maximum range to infinite
maxNSamples = -inf;

if isfield(varargin{1},'format'); bFormat = true; else bFormat = false; end
if isfield(varargin{1},'nBits');  bNBits  = true; else bNBits  = false; end

if bFormat; nFormat = cell(nSignals,1);  end
if bNBits;  nBits   = zeros(nSignals,1); end

% Loop over number of input signals
for ii = 1 : nSignals
    % Check if input is an audio struct 
    if isAudio(varargin{ii})
       % Update audio structure
       varargin{ii} = updateAudio(varargin{ii});
       % Tack max. signal length
       maxNSamples = max(maxNSamples,size(varargin{ii}.data,1));
       % Track sampling rate
       fsIdx(ii) = varargin{ii}.fs;
       % Track number of channels
       nChan(ii) = varargin{ii}.nChannels;
       if bFormat
           nFormat{ii} = varargin{ii}.format;
       end
       if bNBits
           if isfield(varargin{ii},'nBits') && ~isempty(varargin{ii}.nBits)
               nBits(ii) = varargin{ii}.nBits;
           else
               bNBits = false;
           end
       end
    else
        error('Require audio structure as input.')
    end
end

% Check for proper sampling rate 
if any(fsIdx(1) ~= fsIdx)
   error('Sampling rate differs between audio signals!')
else
    fs = fsIdx(1);
end

% Check if number of channels match
if any(nChan(1) ~= nChan)
   error('Number of channels differs between audio signals!')
else
    nChan = varargin{1}.nChannels;
end

% Allocate memory
data = zeros(maxNSamples,nChan);

% Sum all signals
for ii = 1 : nSignals
    currIdx         = 1:size(varargin{ii}.data,1);
    data(currIdx,:) = data(currIdx,:) + varargin{ii}.data(currIdx,:);
end

% Create empty structure
out = struct('label',[],'data',data,'fs',fs);

% Label
if nSignals == 1
    out.label = sprintf('mixture of %i signal',nSignals);
else
    out.label = sprintf('mixture of %i signals',nSignals);
end

% Format
if bFormat
    if ~any(strcmp(nFormat{1},nFormat))
        warning('MATLAB:mixAudio','Format differs between audio signals.')
    else
        out.format = nFormat{1};
    end
end

% Number of bits
if bNBits
    if any(nBits(1) ~= nBits)
        warning('MATLAB:mixAudio',['Number of bits differs between',...
                ' audio signals!'])
    else
        if ~isequal(nBits(1),0)
            out.nBits = nBits(1);
        else
            out.nBits = [];
        end
    end
end

% Refresh fields of audio structure
out = updateAudio(out);


%% ***************************  PLOT MIXTURE  *****************************
% 
% 
% Plot mixture
if nargout == 0
   plotAudio(out);
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