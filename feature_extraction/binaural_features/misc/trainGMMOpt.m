function [gmmFinal,logEM,ubm] = trainGMMOpt(featSpace,classID,C,gmmInit,bVerbose)

%   Developed with Matlab 7.5.0.342 (R2007b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2009-2010 
%              University of Oldenburg and TU/e Eindhoven 
%              tobias.may@uni-oldenburg.de   t.may@tue.nl
%
%   History :
%   v.0.1   2009/11/23
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
narginchk(3, 5);

% Set default values
if nargin < 5 || isempty(bVerbose); bVerbose = false; end

% Determine feature space dimensions
[nObs, nFeatures] = size(featSpace);

% Check consistency
if nObs ~= length(classID)
    error('Mismatch between feature space and the number of class labels');
end

% Find unique number of classes
classIdx = unique(classID);
nClasses = length(classIdx);

% Allocate cell array
logEM = cell(nClasses,1);

% Empty variable
ubm = [];
    
% % Create GMM for each class
% [gmmInit,gmmFinal] = deal(repmat(gmm(nFeatures,C.nClusters,...
%                           C.covarType),[nClasses 1]));
% 
% % Check if UBM should be used
% if strfind(lower(C.method),'ubm')
%     
%     % Copy empty GMM structures to initialize UBM structure
%     ubm = gmmInit;
%     
%     if iscell(C.ubm); szUBM = C.ubm; else szUBM = {C.ubm}; end
%     models = cell(length(szUBM),1);
%     
% 	% Loop over all UBMs
%     for jj = 1 : length(szUBM)
%         if exist(szUBM{jj},'file')
%             % Load all UBMs
%             models{jj} = load(szUBM{jj});
%         else
%             error(['The UBM "',szUBM{jj},'" could not be loaded.']);
%         end
%     end
%     
%     % Initialize UBM parameter
%     paramUBM = configUBM('selectGMMs',C.nClusters);
% else
%     % Empty variable
%     ubm = [];
% end
                      
% Loop over number of classes
for ii = 1 : nClasses
    
    % Find frame indices corresponding to current class
    currClass = classID == classIdx(ii);
    
    % Select feature space of current class
    currFeatSpace = featSpace(currClass,:);
    
    % Report status
    if bVerbose
        fprintf('\nInit. GMM model %i/%i ... ',ii,nClasses)
    end
    
%     % Check if UBM should be used
%     if strfind(lower(C.method),'ubm')
%         llScore = zeros(length(szUBM),1);
%         % Determine best UBM for ii-th class
%         for jj = 1 : length(szUBM)
%             llScore(jj) = sum(log(gmmprob(models{jj}.ubm.GMM,...
%                               currFeatSpace)));
%         end
%         % Determine winning UBM
%         [tmp,idxUBM] = max(llScore); %#ok
%         
%         % Initialize ii-th class 
%         gmmInit(ii) = models{idxUBM}.ubm.GMM;
%         
%         % Clean up
%         clear llScore tmp idxUBM;
%     else
%         % Initialization
%         gmmInit(ii) = initializeGMM(gmmInit(ii),currFeatSpace,C.init);
%     end
            
    % Report status
    if bVerbose
        fprintf('Train GMM model %i/%i ... ',ii,nClasses)
    end

    % Select training method
    switch lower(C.method)
        case 'gmmubm'
            % Adapt a GMM-UBM model for individual classes
            [gmmFinal(ii),logEM{ii},ubm(ii)] = gmmmap(gmmInit(ii),...
                                                      currFeatSpace,...
                                                      C.optEM,...
                                                      paramUBM);
        case 'gmm'
            % Train GMM model for individual classes
            [gmmFinal(ii),logEM{ii}] = gmmem(gmmInit(ii),currFeatSpace,...
                                             C.optEM);
        case 'figueiredo'
            % Automatically select number of Gaussian components
            [gmmFinal(ii),logEM{ii}] = gmmb_fj(gmmInit(ii),currFeatSpace,C);
    end
        
    % Check if the covariance matrix of the ii-th GMM model is badly
    % conditioned. 
    if isequal(lower(C.covarType),'diag')
        if ~checkCovariance(gmmFinal(ii))
            warning(['Covariance matrix of class ',num2str(ii),...
                   ' is badly conditioned.']);
        end
    end
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