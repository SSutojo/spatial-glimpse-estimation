function gmmAz = fitMoG(azEstimate)
%fitMoG   Fit a mixture of Gaussians.
%
%USAGE
%   gmmAz = fitMoG(azEstimate)
%
%INPUT ARGUMENTS
%   azEstimate : azimuth estimate [nFrames x 1]
%
%OUTPUT ARGUMENTS
%   gmmAz : Structure containing mixture of Gaussians (NETLAB format)

%   Developed with Matlab 7.8.0.347 (R2009a). Please send bug reports to:
%
%   Author  :  Tobias May, © 2009
%              TUe Eindhoven and Philips Research
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :
%   v.0.1   2009/08/29
%   v.0.2   2009/11/24 cleaned-up
%   ***********************************************************************

% Check for proper input arguments
error(nargchk(1,1,nargin));


%% **************************  USER PARAMETERS  ***************************

% Determine azimuth range
azMin = -90;
azMax =  90;

% Spread of background model in azimuth
azSpreadBackground = azMax - azMin;

% Define spread of sources according to "gaussianSTD" standard deviations
% of a Gaussian distribution 
% gaussianSTD = 1; % about 68%   data coverage
% gaussianSTD = 2; % about 95%   data coverage
gaussianSTD = 3; % about 99.7% data coverage

% Number of Gaussian components
nClusters = 8;

% Spread of sound source candidates in azimuth (covariance floor)
azSpreadSources = 10;

% Spread of initial sound source candidates in azimuth (should be high so
% that each centre is responsible for a reasonable number of points)
azSpreadInit = 20;

% Merge GMMs if their distance is below 'minDistAz' (shrink)
minDistAz = 10;


% ===================================
% Map spread in azimuth to covariance
% ===================================

% Covariance floor of source candidates
covarFloorSources = ((azSpreadSources/2)/gaussianSTD).^2;

% Initial covariance of GMMs, should be high so that each centre is
% responsible for a reasonable number of points. 
covarInit = ((azSpreadInit/2)/gaussianSTD).^2;

% Mean and covariance floor of background GMM
meanBackground       = 0;
covarFloorBackground = ((azSpreadBackground/2)/gaussianSTD).^2;

% Display training progress
bVerbose = false;

% Plot GMM model adaptation
bPlot = false;

% EM convergence criterion
terminationEM = 1e-5;

% Number of EM iteration steps
nIterEM = 500;

% Number of Kmeans iteration steps
nIterKmeans = 20;

% Delete a Gaussian component if its prior is below this threshold
% thresPrior = 0.05;
thresPrior = 0.0;


%% ***********************  INITIALIZE PARAMETERS  ************************

% EM paramerters
optEM = zeros(1,18);
% Display training progress
optEM(1) = bVerbose;
% EM termination threshold
optEM(3) = terminationEM;

% Covariance flooring
optEM(5) = covarFloorSources;
% Number of iterations for EM training
optEM(14) = nIterEM;

% kMeans
optKMeans = foptions; %#ok (using NETLABs version)
% Number of kMean iterations
optKMeans(14) = nIterKmeans;
% Initialize Kmeans clusters from data
optKMeans(5) = true;


% Create GMM with diagonal covariance matrices including hidden GMM for
% background modeling
gmmAz = gmm(1,nClusters + 1,'diag');

% Check that inputs are consistent
error(consist(gmmAz, 'gmm', azEstimate));


%% ***********************  INITIALIZE GMM MODELS  ************************

% Initialize GMM components (do not change the background model)
gmmIdx = 1:nClusters;

% Initialize GMM component means using kMeans clustering algorithm
[gmmAz.centres(gmmIdx),optKMeans,post] = kmeans(gmmAz.centres(gmmIdx),  ... 
                                                azEstimate,optKMeans); %#ok

% Initialize GMM component variances
for jj = gmmIdx
    % Pick out data points belonging to the jj-th center
    c     = azEstimate(logical(post(:, jj)),:);
    diffs = c - (ones(size(c, 1), 1) * gmmAz.centres(jj, :));
    if size(c, 1) > 0
        gmmAz.covars(jj, :) = sum((diffs.*diffs), 1)/size(c, 1);
        % Replace small entries by covarInit value
        gmmAz.covars(jj, :) = max(gmmAz.covars(jj, :),covarInit);
    else
        % Replace small entries by covarInit value
        gmmAz.covars(jj, :) = covarInit;
    end
end

% Initialize hidden layer
gmmAz.centres(end)   = meanBackground;
gmmAz.covars(end, :) = covarFloorBackground;

bProcess = true;
bMerge   = true;

if bPlot
    h = figure;
    xGrid = azMin:0.5:azMax;
    figure(h);clf;
    plot(xGrid,gmmprob(gmmAz,xGrid.'))
    title(['Mixture of Gaussians consisting of ',...
          num2str(gmmAz.ncentres),' components ',...
          '(including background model)'])
    grid on;
    xlim([azMin azMax])
    xlabel('Azimuth (deg)')
    ylabel('Sound source evidence')
    pause(0.5)
end


%% **********************  PERFORM GMM ADAPTATION  ************************

while bProcess
    
    if bMerge
        bProcess = true;
    else
        bProcess = false;
    end
    
    % Adapt GMM structure
    gmmAz = adaptGMM(gmmAz,azEstimate,optEM,covarFloorBackground);
    
    % Plot GMM model
    if bPlot
        figure(h);clf;
        plot(xGrid,gmmprob(gmmAz,xGrid.'))
        title(['Mixture of Gaussians consisting of ',...
              num2str(gmmAz.ncentres),' components ',...
               '(including background model)'])
        grid on;
        xlim([azMin azMax])
        xlabel('Azimuth (deg)')
        ylabel('Sound source evidence')
        drawnow;
    end
    
    % Sort centres
    [centresRanked,sortIdx] = sort(gmmAz.centres(1:end-1),'descend');
    
    % Check source distance
    distMoG = abs(diff(centresRanked));
    
    % Check if GMM components should be merged
    mergeIdx = distMoG < minDistAz;
    
    % Merge GMM components which are closest together ...
    if any(mergeIdx)
        
        idx2Merge = find(mergeIdx);
        [m,mIdx]  = min(distMoG(idx2Merge));
        mergeIdx  = idx2Merge(mIdx);
       
        m1 = sortIdx(mergeIdx);
        m2 = sortIdx(mergeIdx+1);
           
        priorM  =  gmmAz.priors(m1) + gmmAz.priors(m2);        
        centreM = (gmmAz.priors(m1) * gmmAz.centres(m1) + ...
                   gmmAz.priors(m2) * gmmAz.centres(m2)) / priorM;
        covarM  = (gmmAz.priors(m1) * gmmAz.covars(m1)  + ...
                   gmmAz.priors(m2) * gmmAz.covars(m2)) / priorM;
        
        % Shrink number of components
        gmmAz.ncentres = gmmAz.ncentres - 1;
        gmmAz.nwts     = gmmAz.ncentres + gmmAz.ncentres*gmmAz.nin +...
                         gmmAz.ncentres * gmmAz.nin;
         
        rmIdx    = max(m1,m2);
        storeIdx = min(m1,m2);
        
        gmmAz.centres(rmIdx) = [];
        gmmAz.covars(rmIdx)  = [];
        gmmAz.priors(rmIdx)  = [];
        
        gmmAz.centres(storeIdx) = centreM;
        gmmAz.covars(storeIdx)  = covarM;
        gmmAz.priors(storeIdx)  = priorM;
        
        % Rescale priors
        % gmmAz.priors = ones(1,gmmAz.ncentres) ./ gmmAz.ncentres;
        
        bMerge = true;
    else
        % Deactivate merging
        bMerge = false;
    end
end

% Remove Gaussian components which are below this threshold
rmIdx = gmmAz.priors(1:end-1) < thresPrior;

if any(rmIdx)
    % Shrink number of components
    gmmAz.ncentres = gmmAz.ncentres - sum(rmIdx);
    gmmAz.nwts     = gmmAz.ncentres + gmmAz.ncentres*gmmAz.nin +...
                     gmmAz.ncentres * gmmAz.nin;
    
    gmmAz.centres(rmIdx) = [];
    gmmAz.covars(rmIdx)  = [];
    gmmAz.priors(rmIdx)  = [];
    
    % Re-Adapt GMM structure
    gmmAz = adaptGMM(gmmAz,azEstimate,optEM,covarFloorBackground);
end


%% ***************************  SHOW RESULT  ******************************

if  bPlot || nargout == 0
    % close(h);
    
    azGridHist = azMin:2.5:azMax;
    azHist     = histc(azEstimate(:),azGridHist);
    azHist     = azHist ./ numel(azEstimate);

    figure;hold on;
    h(1)  = bar(azGridHist,azHist,1);
    xGrid = azMin:0.5:azMax;
    h(2)  = plot(xGrid,gmmprob(gmmAz,xGrid.'),'Linewidth',2,'Color','r');
    legend(h,{'azimuth histogram' 'MoG'})
    title(['Mixture of Gaussians consisting of ',...
           num2str(gmmAz.ncentres),' components ',...
           '(including background model)'])
    grid on;
    xlim([azMin azMax])
    xlabel('Azimuth (deg)')
    ylabel('Sound source evidence')
end

% =========================================================================
% SUBFUNCTION
% =========================================================================
% 
function gmmAz = adaptGMM(gmmAz,azEstimate,opt,covarFloorBackground)
%adaptGMM   Adapt GMM parameter.

% Get dimension of input data
[ndata, xdim] = size(azEstimate); %#ok

% Sort out the options
if (opt(14))
    niters = opt(14);
else
    niters = 100;
end

display = opt(1);
store   = 0;
if (nargout > 2)
    store = 1;	% Store the error values to return them
    errlog = zeros(1, niters);
end
test = 0;
if opt(3) > 0.0
    test = 1;  % Test log likelihood for termination
end

if opt(5) > 0
    if display >= 0
        disp('check_covars is on');
    end
    check_covars = 1;	   % Ensure that covariances don't collapse
    MIN_COVAR    = opt(5); % Minimum singular value of covariance matrix
else
    check_covars = 0;
end

% Main loop of algorithm
for n = 1:niters
    
    % Calculate posteriors based on old parameters
    [post, act] = gmmpost(gmmAz, azEstimate);
    
    % Calculate error value if needed
    if (display || store || test)
        prob = act*(gmmAz.priors)';
        % Error value is negative log likelihood of data
        e = - sum(log(prob));
        if store
            errlog(n) = e;
        end
        if display > 0
            fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
        end
        if test
            if (n > 1 && abs(e - eold) < opt(3))
                opt(8) = e; %#ok
                return;
            else
                eold = e;
            end
        end
    end
    
    % Adjust the new estimates for the parameters
    new_pr = sum(post, 1);
    new_c = post' * azEstimate;
    
    % Now move new estimates to old parameter vectors (including background
    % GMM) 
    gmmAz.priors = new_pr ./ ndata;
    
    % Do not update means of hidden layer
    gmmAz.centres(1:end-1) = new_c(1:end-1) ./ (new_pr(1:end-1)' * ...
                             ones(1, gmmAz.nin));
    
    % Update covariance
    for j = 1:gmmAz.ncentres
        diffs = azEstimate - (ones(ndata, 1) * gmmAz.centres(j,:));
        if j == gmmAz.ncentres
            % Hidden layer
            gmmAz.covars(j,:) = sum((diffs.*diffs).*(post(:,j)*ones(1,...
                                gmmAz.nin)), 1)./new_pr(j);
            % Apply flooring
            gmmAz.covars(j,:) = max(covarFloorBackground,gmmAz.covars(j,:));
        else
            % Source candidates
            gmmAz.covars(j,:) = sum((diffs.*diffs).*(post(:,j)*ones(1,...
                                    gmmAz.nin)), 1)./new_pr(j);
        end
    end
    if check_covars
        % Ensure that no covariance is too small
        for j = 1:gmmAz.ncentres
            gmmAz.covars(j,:) = max(gmmAz.covars(j,:),MIN_COVAR);
        end
    end
end