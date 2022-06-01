function predicted = nntest(nn, x)

% Ensemble size
nNN = numel(nn);

% Set training flag
bTrain = false;

if nNN > 1
       
    % Allocate memory
    predicted = zeros(size(x,1),size(nn(1).W{end},1));
    
    % Loop over the enseble
    for ii = 1 : nNN
        
        % Forward pass
        A = nnff(nn(ii), x, bTrain);
        
        % Return activation of output layer
        predicted = predicted + A.states.A{end};
    end
    
    % Normalization
    predicted = predicted / numel(nn);
    
else
    
    % Forward pass
    nn = nnff(nn, x, bTrain);
    
    % Return activation of output layer
    predicted = nn.states.A{end};
    
end