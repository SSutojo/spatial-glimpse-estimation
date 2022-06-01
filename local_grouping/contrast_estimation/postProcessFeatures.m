function [FEATURES,LABELS,M] = postProcessFeatures(FEATURES,M,fileIdx,LABELS)


%% CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Identify modus
if nargin == 4 
    if isfield(M,'mapping') && ~isempty(M.mapping)
        mode = 'valid';
    else
        mode = 'train';
    end
elseif nargin == 2
    mode    = 'test';
    fileIdx = [1 size(FEATURES,1)];
else
    error('Wrong usage.')
end


%% LIMIT TRAINING DATA
% 
% 
% Total number of files
nFiles = size(fileIdx,1);

% Select subset of training data
if strcmpi(mode,'train') && M.opt.post.trainRatio < 1
    nFiles = round(M.opt.post.trainRatio * nFiles);
    LABELS = LABELS(1:fileIdx(nFiles,2),:);
end


%% COMPRESSION
% 
% 
% Apply compression and accumulate features
FEATURES = compress(FEATURES(1:fileIdx(nFiles,2),:),M.opt.post.compress{:});

            
%% EXPAND TEMPORAL FEATURE CONTEXT
% 
% 
% Expand temporal context
if any(M.opt.post.context.context > 0)
    
    % Numer of frames [past future]
    K = M.opt.post.context.context;
    
    % New feature dimension
    FNew = numel(-K(1):K(1));
    
    % Allocate memory
    FEAT_TMP = zeros(size(FEATURES,1), FNew * size(FEATURES,2));
    
    % File-based processing
    for ii = 1 : nFiles
        idx = (fileIdx(ii,1)):fileIdx(ii,2);
        
        % Loop over feature dimension
        for ff = 1 : size(FEATURES,2)
            
            % New feature indices
            idxF = (1:FNew) + (ff-1) * (FNew + M.opt.post.context.bAppend);
            
            % Expand temporal context
            FEAT_TMP(idx,idxF) = integrateContext(FEATURES(idx,ff).',...
                K,M.opt.post.context.bDCT,...
                M.opt.post.context.orderDCT + ...
                M.opt.post.context.bRemoveDC, ...
                M.opt.post.context.bRemoveDC).';
            
            if M.opt.post.context.bAppend
                % Append feature refelcting temporal context 
                FEAT_TMP(idx,idxF(end) + 1) = FEATURES(idx,ff);
            end
        end
    end
    
    FEATURES = FEAT_TMP;
    clear FEAT_TMP;
end


%% FEATURE NORMALIZATION
% 
% 
% Supported normalization strategies
strNorm = {'mean' 'var' 'meanvar' 'max' 'range' 'heq'};

% Local normalization (file-based)
if any(strcmpi(M.opt.post.normalize.local.normMethod,strNorm))
    for ii = 1 : nFiles
        idx = (fileIdx(ii,1)):fileIdx(ii,2);
        FEATURES(idx,:) = normalizeData_dnn(FEATURES(idx,:),...
            M.opt.post.normalize.local.normMethod,1);
    end
end

% Global normalization
if any(strcmpi(M.opt.post.normalize.global.normMethod,strNorm))
    if strcmpi(mode,'train')
        [FEATURES,M.mapping] = normalizeData_dnn(FEATURES,...
            M.opt.post.normalize.global.normMethod,1);
    else
        FEATURES = scaleData(FEATURES,M.mapping);
    end
else
    M.mapping = false;
end
