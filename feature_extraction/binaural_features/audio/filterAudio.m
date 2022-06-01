function [in,fObj] = filterAudio(in,fObj)
%filterAudio   Perform digital filtering of data/audio objects.
%
%USAGE 
%   [in,fObj] = filterAudio(in,fObj)
%
%INPUT ARGUMENTS
%     in : either audio object (see genAudio) or data matrix arranged as 
%          [nSamples x nChannels]
%   fObj : filter object(s)
%
%OUTPUT ARGUMENTS
%     in : filtered audio object or data matrix
%   fObj : filter object with updated filter states
%
%NOTE
%   If multiple filter objects are cascaded in 'parallel' or in 'series',
%   the dimension of the filtered data will be equal to the input data. If
%   multiple filter objects should be processed individually, a third
%   dimension will be added which corresponds to the number of filter
%   objects. See "cascadeFilter" for more details.
% 
%EXAMPLE
% 
%   See also genFilterObj, plotFilter, calcFilterResponse and 
%            cascadeFilter.


%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2007-2009 
%              TUe Eindhoven 
%              t.may@tue.nl   
%
%   History :  
%   v.0.1   2007/08/05
%   v.0.2   2007/08/21 added option to cascade of filter objects 
%   v.0.3   2008/09/15 added audio structure support
%   v.0.4   2009/10/20 cleaned up
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
narginchk(2, 2);

% Unpack and check filter objects 
fObj = unpackFilterObj(fObj);

% Number of filter objects
nFilter = length(fObj);

% Check if "in" is an audio structure
if isAudio(in)
   fs      = in.fs;
   data    = in.data;
   bAudio  = true;
else
    % Get dimension of input signal
    data   = in;
    bAudio = false;
end

% Check dimensionality of data
[nSamples,nChannels,hdim] = size(data);

% Check if data is 2 dimensional
if hdim > 1
    error('"IN" must be either 1 or 2 dimensional.')
end

% Multiple filter objects
if nFilter > 1
    % Check for consistency
    if ~all([fObj.fs] == fObj(1).fs)
        error('Sampling frequency is not consistent across filter objects');
    end
    if ~all(strcmp(fObj(1).cascade,{fObj.cascade}))
       error('"Cascade" option is not consistent across filter objects');
    end
    if ~all(strcmp(fObj(1).structure,{fObj.structure}))
        error('Filter structure is not consistent across filter objects');
    end
end

% Check for consistent sampling frequency
if bAudio && ~isequal(fs,fObj(1).fs)
    error('Sampling frequency mismatch between audio and filter object.')
end

% Allocate memory
switch lower(fObj(1).cascade)
    case {'parallel' 'series'}
        out = zeros(nSamples,nChannels);
    case {''  'individual'}
        out = zeros(nSamples,nChannels,nFilter);
    otherwise
        error(['Cascade option "',lower(fObj(1).cascade),...
               '" is not recognized.'])
end

% Select filter structure
switch fObj(1).structure
    % case 'Direct-Form II'
    %     szFilter = 'filterDF2';

    % TODO % Add filter structure Direct-Form II

    case 'Direct-Form II Transposed'
        szFilter = 'filter';
    otherwise
    error(['Filter structure "',fObj(1).structure,'" is not recognized.'])
end


%% *****************************  FILTERING  ******************************
% 
% 
% Loop over number of filter objects
for ii = 1 : nFilter
    
    % *******************  CHECK FILTER COEFFICIENTS  *********************
    
    % Get dimension of filter coefficients
    dimB = size(fObj(ii).b);
    dimA = size(fObj(ii).a);

    if dimB(1) == 1 && dimB(2) > 1; fObj(ii).b = fObj(ii).b(:);
        % Get dimension of filter coefficients
        dimB = size(fObj(ii).b);
    end
    if dimA(1) == 1 && dimA(2) > 1; fObj(ii).a = fObj(ii).a(:);
        % Get dimension of filter coefficients
        dimA = size(fObj(ii).a);
    end

    % **********************  CHECK FILTER STATES  ************************
    
    % Check if filter states are initialized
    if isempty(fObj(ii).states) 
        % Initialize filter states
        fObj(ii).states = zeros(max(length(fObj(ii).b),...
                                length(fObj(ii).a))-1,nChannels);
    else
        if max(dimB(1),dimA(1)) - 1 ~= size(fObj(ii).states,1)
            error(['Dimension mismatch between the filter ',...
                   'coefficients and the filter states.']);
        end
    end
    
    % Extend filter states to the number of audio channls
    if size(fObj(ii).states,2) ~= nChannels
        % Check if filter states are zero
        if sum(fObj(ii).states) == 0 
            % Replicate filter states
            fObj(ii).states = repmat(fObj(ii).states,[1 nChannels]);
        else
            error(['Dimension mismatch between the filter states and ' ,...
                   'the number of audio channels! The Initialized '    ,...
                   'filter states are non-zero and are therefore '     ,...
                   'assumed to be data-dependent. However, the second ',...
                   'dimension of the filter states does not '          ,...
                   'correspond to the number of audio channels.'])  
        end
    end

    
    % ***********************  PERFORM FILTERING  *************************
    
    % Perform filtering
    switch fObj(ii).cascade
        case 'parallel'
            % Sum the output of all filter objects
            [outP,fObj(ii).states] = feval(szFilter,fObj(ii).b,...
                                           fObj(ii).a,data,...
                                           fObj(ii).states);
            out = out + outP;
        case 'series'
            % Cascade all filter objects in series
            [out,fObj(ii).states] = feval(szFilter,fObj(ii).b,...
                                          fObj(ii).a,data,fObj(ii).states);
            data = out;
        case {'' 'individual'}
            % Create output for each filter object
            [out(:,:,ii),fObj(ii).states] = feval(szFilter,fObj(ii).b,...
                                                  fObj(ii).a,data,...
                                                  fObj(ii).states);
    end
end


%% **********************  RETURN FILTERED SIGNAL  ************************
% 
% 
% Create output
if bAudio
    % Copy signal
    in.data = out;
    % Update audio structure
    in = updateAudio(in);
else
    in = out;
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