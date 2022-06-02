% DEMGMMC Demo of GMM used for classification: version to produce figure

% Work with 2-d data
nin = 2;
% Fix the seeds
rand('state', 4231);
randn('state', 4231);

% 
% Generate the data: this follows demmlp2
% 
n=200;

% Set up mixture model: 2d data with three centres
% Class 1 is first centre, class 2 from the other two
mix = gmm(nin, 3, 'full');
mix.priors = [0.5 0.25 0.25];
mix.centres = [0 -0.1; 1 1; 1 -1];
mix.covars(:,:,1) = [0.625 -0.2165; -0.2165 0.875];
mix.covars(:,:,2) = [0.2241 -0.1368; -0.1368 0.9759];
mix.covars(:,:,3) = [0.2375 0.1516; 0.1516 0.4125];

[data, label] = gmmsamp(mix, n);

% Now fit a mixture of Gaussians to each class. Use spherical 
% Gaussians as we assume we don't know the data structure.
ncentres = 7;
mix1 = gmm(nin, ncentres, 'spherical'); 
mix2 = gmm(nin, ncentres, 'spherical');

% Initialise models
options = foptions;
options(14) = 5;	% Just use 5 iterations of k-means in initialisation
% Initialise the model parameters from the data
mix1 = gmminit(mix1, data(label==1, :), options);
mix2 = gmminit(mix2, data(label>1, :), options);

% Now train each model
options(1) = 1;  % Print error values
options(5) = 1;  % Prevent collapse of variances
options(14) = 15; % Number of iterations
mix1 = gmmem(mix1, data(label==1, :), options);
mix2 = gmmem(mix2, data(label>1, :), options);

disp('Press any key to continue');
pause

% Generate test data
[test, test_label] = gmmsamp(mix, n);
% Priors are equal and normalising constant is the same for both 
% classes, so cannot affect classification
test_probs = [gmmprob(mix1, test) gmmprob(mix2, test)];
% Convert to 1 of N encoding
target = [test_label==1 test_label>1];

fh1 = conffig(test_probs, target);

% Now produce some interesting plots
MarkerSize = 4;
fh2 = figure;
% Extract different classes from test data
test_c1 = test(test_label==1, :);
test_c2 = test(test_label>1, :);
p1 = plot(test_c1(:, 1), test_c1(:, 2), 'ko', ...
  'MarkerSize', MarkerSize);
hold on;
plot(test_c2(:, 1), test_c2(:, 2), 'kx', 'MarkerSize', MarkerSize);
MarkerSize=10;
plot(mix1.centres(:, 1), mix1.centres(:, 2), 'ks', ...
    'MarkerSize', MarkerSize, 'MarkerFaceColor', 'k');
plot(mix2.centres(:, 1), mix2.centres(:, 2), 'kv', ...
    'MarkerSize', MarkerSize, 'MarkerFaceColor', 'k');
%legend('Class 1', 'Class 2', 'Mix 1', 'Mix2');
hold on;

% Now plot decision boundaries
x0 = min(test(:,1));
x1 = max(test(:,1));
y0 = min(test(:,2));
y1 = max(test(:,2));
dx = x1-x0;
dy = y1-y0;
expand = 5/100;			% Add on 5 percent each way
x0 = x0 - dx*expand;
x1 = x1 + dx*expand;
y0 = y0 - dy*expand;
y1 = y1 + dy*expand;
resolution = 100;
step = dx/resolution;
xrange = [x0:step:x1];
yrange = [y0:step:y1];
% 					
% Generate the grid
% 
[X Y]=meshgrid([x0:step:x1],[y0:step:y1]);
% Posterior of mix gives optimal decision boundary
post = gmmpost(mix, [X(:) Y(:)]);
p1_x = reshape(post(:, 1), size(X));
p2_x = reshape(post(:, 2) + post(:, 3), size(X));
[c, h1] = contour(xrange,yrange,p1_x,[0.5 0.5],'k-');
axis([x0 x1 y0 y1])
% Now compute our decision boundary; it's where the two probabilities
% probabilities from mix1 and mix2 are equal (again, assuming equal
% priors
prob1 = gmmprob(mix1, [X(:) Y(:)]);
prob2 = gmmprob(mix2, [X(:) Y(:)]);
prob_diff = reshape(prob1-prob2, size(X));
[c, h] = contour(xrange, yrange, prob_diff, [0 0], 'k-.');
set(h, 'LineWidth', 2);
print -deps demgmmc.eps;

% Finally, compute the confusion matrix for the optimal rule
postrule = gmmpost(mix, test);
fh3 = conffig([postrule(:, 1) (postrule(:, 2)+postrule(:, 3))], ...
  target);
