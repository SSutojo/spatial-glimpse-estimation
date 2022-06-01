function [fSpace,idx] = stackFrames(fSpaceIn,cDim,bTrain,labelsIn)
%stackFrames   Stack matrix across neighboring frames.
% 
%USAGE
%            fSpace = stackFrames(fSpaceIn,cDim)
%            fSpace = stackFrames(fSpaceIn,cDim,bTrain,labelsIn)
%   [fSpace,labels] = stackFrames(fSpaceIn,cDim,bTrain,labelsIn)
% 
%INPUT ARGUMENTS
%     fSpaceIn : feature matrix [fDim x nFrames]
%         cDim : number of neighboring frames across which the input matrix
%                should be stacked. The first dimension refers to the
%                number of past frames, whereas the second dimension
%                defines the number of future frames (default, cDim = 5)
%       bTrain : select training or testing mode, if bTrain == true
%                (training mode) the output matrix is trimmed so that only
%                frames with the complete context are retained, if bTrain
%                == false (testing mode), the first and the last feature
%                frames are repeated to fill the temporal context at the
%                beginning and the end (default, bTrain = true)
%     labelsIn : 2D label matrix which is trimmed according to training or
%                testing mode [lDim x nFrames]
% 
%OUTPUT ARGUMENTS
%       fSpace : 2D feature space [fDim * cDim x nFramesTrimmed]
%       labels : 2D label matrix [lDim x nFramesTrimmed]

%   Developed with Matlab 9.1.0.441655 (R2016b). Please send bug reports to
%   
%   Author  :  Tobias May, © 2016
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2016/11/03
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 2 || isempty(cDim);   cDim   = 5;    end
if nargin < 3 || isempty(bTrain); bTrain = true; end

% Check dimContext
if any(~rem(cDim,1)==0) || any(cDim < 0)
    error('"dimContext" must be nonzero integers.')
end


%% STACK FEATURE SPACE ACROSS NEIGHBORING FRAMES
%
%
% Get dimensions
[fDim,nFrames] = size(fSpaceIn);

% Generate frame indices
switch length(cDim)
    case 1
        rowIdx = (-cDim(1):0)';
    case 2
        rowIdx = (-cDim(1):cDim(2))';
    otherwise
        error('Context must have either 1 or 2 dimensions.')
end

% Find first valid frame index 
startIdx = find(rowIdx==1)-1;

% Set default value (distance to first causal frame)
if isempty(startIdx)
    startIdx = argmax(rowIdx);
%     startIdx = argmax(rowIdx) + 1;
end

% Dimension of across-frame context
nFramesContext = numel(rowIdx);

% Create frame indices
colIdx   = 1 + (0:(nFrames-1));
frameIdx = rowIdx(:,ones(1,nFrames)) + colIdx(ones(numel(rowIdx),1),:);

if bTrain
    % Limit range
    frameIdx = frameIdx(:,startIdx:end-(startIdx-1));

    % Allocate memory
    fSpace = zeros(fDim * nFramesContext, nFrames - 2 * (startIdx - 1));
else
    % Limit range
    frameIdx = min(max(frameIdx,1),nFrames);

    % Allocate memory
    fSpace = zeros(fDim * numel(rowIdx), nFrames);
end

% Loop over feature space dimension
for kk = 1 : fDim

    % Get kk-th feature dimension
    fSpaceTmp = fSpaceIn(kk,:).';
    
    % Indices of stacked feature space
    fSpaceIdx = (1:nFramesContext) + (kk-1) * nFramesContext;
    
    % Stack context
    fSpace(fSpaceIdx,:) = fSpaceTmp(frameIdx);
end

idx = startIdx:nFrames-(startIdx-1);

% Select proper labels
if exist('labelsIn','var')
    labels = labelsIn(:,startIdx:end-(startIdx-1));
end
