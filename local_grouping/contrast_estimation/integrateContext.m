function y = integrateContext(x,cDim,bDCT,N,bRemoveDC)

if nargin < 5 || isempty(bRemoveDC); bRemoveDC = false; end
if nargin < 4 || isempty(N);         N         = 5;     end
if nargin < 3 || isempty(bDCT);      bDCT      = true;  end

% Stack context
y = stackFrames(x,cDim,false);

if bDCT
    D = size(y,1);
    
    % DCT matrix
    mDCT = cos(pi * (0:N-1)' * ((1:D) - 0.5) / D);
    
    % Apply DCT
    y = mDCT * y;
    
    if bRemoveDC
        y = y(2:end,:);
    end
end