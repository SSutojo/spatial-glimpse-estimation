function out = updateAudio(in,label)
%updateAudio   Updates structure fields of an audio object.
%   An audio object is a MATLAB structure with at least two fiels,namely
%   "DATA" and "FS". 
% 
%USAGE
%   OUT = updateAudio(IN)
%   OUT = updateAudio(IN,LABEL)
%   
%INPUT ARGUMENTS
%      IN : audio object 
%   LABEL : label describing audio signal (default, LABEL = [])
% 
%OUTPUT ARGUMENTS
%     OUT : updated audio structure 
% 
%   This function updates and adds the following fields: 
%   out. 
%       label     - label (e.g. wave filename or signal name)
%       nBits     - number of bits
%       nSamples  - length of audio signal in samples
%       nChannels - the number of audio channels
%       range     - start and end point in samples [start end]
%       format    - file format 
% 
%NOTE
%   If the number of bits is not specified, the field will be removed from
%   the audio object.
% 
%EXAMPLE
%   % Load audio material
%   load chirp
%   % Create audio object
%   in = struct('data',y,'fs',Fs,'label','chirp');
%   % Update audio object
%   updateAudio(in)
% 
%   See also genAudio, isAudio and plotAudio.

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2008-2010
%              TUe Eindhoven 
%              t.may@tue.nl     
% 
%   History :  
%   v.0.1   2008/05/15
%   v.0.2   2009/10/20 cleaned up
%   ***********************************************************************

% Check for proper input arguments
narginchk(1, 2);

% Set default value
if nargin < 2 || isempty(label); label = ''; 
    overWriteLabel = false;
else
    overWriteLabel = true;
end

% Check if IN is an audio structure
if ~isAudio(in);
   error('"IN" must be an audio structure.') 
end

% Create empty audio structure
out = struct('label',[],'data',[],'fs',[],'nBits',[],'nSamples',[],...
             'nChannels',[],'range',[],'format',[]);

% Determine audio dimensions
dim = size(in.data);
       
% Update audio structure information 
out.nSamples  = dim(1);
out.nChannels = dim(2);

% Copy DATA and FS
out.data = in.data;
out.fs   = in.fs;

% Update audio range if field is empty or range does not match with number
% of samples
if ~isfield(in,'range') || isempty(in.range) || diff(in.range)+1 ~= dim(1)
    out.range = [1 dim(1)] ;
else
    out.range = in.range;
end

% Label information
if overWriteLabel
    out.label = label;
else
    if isfield(in,'label')
        % Copy field information
        out.label = in.label;
    else
        out.label = label;
    end    
end

% Estimate number of Bits
if (~isfield(in,'nBits') || isempty(in.nBits)) 
    
    % Remove field
    out = rmfield(out,'nBits');
    
%     % Guess data format
%     switch lower(class(in.data))
%         case {'double' 'int64' 'uint64'}
%             out.nBits = 64;
%         case {'single' 'int32' 'uint32'}
%             out.nBits = 32;
%         case {'int16' 'uint16'}
%             out.nBits = 16;            
%         case {'uint8' 'int8'}
%             out.nBits = 8;
%         otherwise
%             error(['Data format "',lower(class(in.data)),...
%                    '" is not recognized.'])
%     end
else
    % Copy field information
    out.nBits = in.nBits;
end

% Set format 
if ~isfield(in,'format') || isempty(in.format)
    % Empty format 
    out.format = '';
else
    % Copy field information
    out.format = in.format;
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