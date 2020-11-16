% standard coefficient solver
%implementation of the following : https://web.media.mit.edu/~crtaylor/calculator.html
stencilPoints = [-4,-3,-2,-1,0,1,2];
desiredDerivitive = 1; %1 is df/dx
bVector = zeros(length(stencilPoints),1);
bVector(desiredDerivitive+1) = factorial(desiredDerivitive);
AMatrix = zeros(length(stencilPoints),length(stencilPoints));
for i = 1: length(stencilPoints)
    
    AMatrix(i,:) = stencilPoints.^(i-1);
end 

AMatrix
bVector
coefficients = AMatrix\bVector
