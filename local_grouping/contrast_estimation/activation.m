function y = activation(x,type,bDerivative)
%activation   Set of activation functions used for neural networks
%   The forward pass computes the activation, while the backward pass
%   requires the derivative of the activation function. Since the
%   activations during the forward pass are stored, the calculation of the
%   derivative assumes that 


if nargin < 3 || isempty(bDerivative); bDerivative = false; end


%% APPLY ACTIVATION FUNCTION
% 
% 
% Select type of activation
switch lower(type)
    
    case 'linear'
        % Linear unit
        if bDerivative
            y = ones(size(x));
        else
            y = x;
        end
        
     case 'sigmoid'
        % Sigmoid
        if bDerivative
            y = x .* (1 - x);
        else
            y = 1 ./ (1 + exp(-x));
        end
        
    case 'tanh'
        if bDerivative
            y = 1 - x.^2;
        else
            y = tanh(x);
        end
        
    case 'tanh_opt'
        if bDerivative
            y = 1.7159 * 2/3 * (1 - 1/(1.7159)^2 * x.^2);
        else
            y = 1.7159 * tanh(2/3 .* x);
        end
        
    case 'elu'
        % Exponential linear unit (ELU) with alpha = 1
        if bDerivative
            bBelowZero = x < 0;
            y = ones(size(x));
            y(bBelowZero) = x(bBelowZero) + 1;
        else
            bBelowZero = x < 0;
            y = x;
            y(bBelowZero) = exp(x(bBelowZero)) - 1;
        end        
        
    case 'relu'
        % Rectified linear unit (ReLU)
        if bDerivative
            y = double(x > 0);
        else
            y = max(x,0.0);
        end
        
    case 'relu6'
        % Rectified linear unit 6 (ReLU6)
        if bDerivative
            bWithin = x > 0 & x < 6;
            y = zeros(size(x));
            y(bWithin) = 1;
        else
            y = min(max(x,0.0),6);
        end
        
    case 'lrelu'
        % Leaky rectified linear unit (LReLU) with alpha = 0.01
        if bDerivative
            y = 0.01 * ones(size(x));
            y(x > 0) = 1;
        else
            y = max(x,0) + 0.01 .* min(x,0);
        end
        
    case 'swish'
        % Swish activation function
        if bDerivative
            y = x + activation(x,'sigmoid',false) .* (1 - x);
        else
            y = x .* (1 + exp(-x)).^-1;
        end
        
    case 'modu'
        % Modulus (ModU) activation function
        if bDerivative
            y = sign(x);
        else
            y = abs(x);
        end
        
    case 'cube'
        % Cube activation [1]
        if bDerivative
            y = 2 .* x.^2;
        else
            y = x.^3;
        end
        
    case 'softmax'
        if bDerivative
            error('Derivative of "softmax" is currently not supported!')
        else
            y = exp(bsxfun(@minus, x, max(x,[],2)));
            y = bsxfun(@rdivide, y, sum(y, 2));
        end
        
    otherwise
        error('Activation ''%s'' is not supported.',lower(type))
end



%% PLOT ACTIVATION FUNCTION
% 
% 
% If no output is specified
if nargout == 0
    
    % Sort data to get a proper input/output curve
    [x2plot,sortIdx] = sort(x,'ascend');
    y2plot = y(sortIdx);
    
    % Plot input/output curve
    figure;
    plot(x2plot,y2plot,'linewidth',1.5,'color',[0 0.3895 0.9712])
    xlabel('Input')
    ylabel('Output')
    grid on;
end
        