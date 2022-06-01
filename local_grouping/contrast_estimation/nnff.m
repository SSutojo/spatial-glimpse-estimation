function nn = nnff(nn, x, bTrain)
%NNFF performs a feedforward pass
% nn = nnff(nn, x) returns an neural network structure with updated
% layer activations

N = nn.architecture.nLayers;
dimB = size(x, 1);

% Input layer + bias
nn.states.A{1} = [ones(dimB,1) x];

% Set derivative flag for activation functions
bDerivative = false;

% Feedforward pass through the hidden layers
for ii = 2 : N-1
    
    % Forward pass A = X * W + b
    nn.states.A{ii} = nn.states.A{ii - 1} * nn.W{ii - 1}';
    
    % If no batch norm or batch norm after activation should be used
    if ~nn.regularization.bBatchNorm || (nn.regularization.bBatchNorm ...
            && ~nn.regularization.bBNBeforeActivation)
        
        % Activation function    
        nn.states.A{ii} = activation(nn.states.A{ii},...
            nn.architecture.activationHidden,bDerivative);
    end
    
    % Batch norm layer
    if nn.regularization.bBatchNorm
        if bTrain
            % Mini-batch mean
            batchMean = mean(nn.states.A{ii},1);
            
            % Centered batch
            xc = (nn.states.A{ii} - repmat(batchMean, dimB, 1));
            
            % Mini-batch variance and standard deviation
            batchVar = mean(xc.^2, 1);
            invstd = 1./sqrt(batchVar + nn.regularization.epsilonBN);
            
            % Normalized batch
            xn = xc .* repmat(invstd, dimB, 1);
            
            % Scale and shift normalized batch
            nn.states.A{ii} = repmat(nn.gamma{ii}, dimB, 1) .* xn + ...
                repmat(nn.beta{ii}, dimB, 1);
            
            % Update mean and variance using exponential averaging
            nn.mu{ii} = nn.regularization.alphaBN * nn.mu{ii} + ...
                (1-nn.regularization.alphaBN) * batchMean;
            nn.var{ii} = nn.regularization.alphaBN * nn.var{ii} + ...
                (1-nn.regularization.alphaBN) * batchVar;
            
            % Store parameters for backpropagation
            nn.states.batchNorm.xn{ii} = xn;
            nn.states.batchNorm.invstd{ii} = invstd;
        else
            % Normalize batch
            xn = (nn.states.A{ii} - repmat(nn.mu{ii}, dimB, 1)) ...
                ./ repmat(sqrt(nn.var{ii} + ...
                nn.regularization.epsilonBN), dimB, 1);
            
            % Scale and shift normalized batch
            nn.states.A{ii} = repmat(nn.gamma{ii}, dimB, 1) ...
                .* xn + repmat(nn.beta{ii}, dimB, 1);
        end
    end
    
    % Only if batch norm should be used before the activation
    if nn.regularization.bBatchNorm ...
            && nn.regularization.bBNBeforeActivation
        
        % Activation function
        nn.states.A{ii} = activation(nn.states.A{ii},...
            nn.architecture.activationHidden,bDerivative);
    end
    
    % Dropout
    if (nn.regularization.dropoutHidden > 0) && ...
            (~nn.regularization.bDropoutLast || ...
            (nn.regularization.bDropoutLast && ii == N-1))
        if bTrain
            % Create and store dropout mask
            nn.states.mask{ii} = (rand(size(nn.states.A{ii})) > ...
                nn.regularization.dropoutHidden);
            
            % Inverted dropout: Scale mask during training, such that no
            % scaling is required during testing
            if nn.regularization.bInvertedDropout
                nn.states.mask{ii} = nn.states.mask{ii} / ...
                    (1 - nn.regularization.dropoutHidden);
            end
            
            % Apply mask
            nn.states.A{ii} = nn.states.A{ii} .* nn.states.mask{ii};
        else
            if ~nn.regularization.bInvertedDropout
                % Scale activation during testing
                nn.states.A{ii} = nn.states.A{ii} .* ...
                    (1 - nn.regularization.dropoutHidden);
            end            
        end
    end
    
    % Add the bias term
    nn.states.A{ii} = [ones(dimB,1) nn.states.A{ii}];
end

% Output layer
nn.states.A{N} = activation(nn.states.A{N - 1} * nn.W{N - 1}',...
    nn.architecture.activationOut,bDerivative);

