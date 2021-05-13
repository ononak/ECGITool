function [ x_reg, errNorm, constNorm ] = Tikhonov( A, R, y, lambda)
% TIKHONOV 
%  Computes the Tikhonov regularized solution
%  min { || A x - y ||^2 + lambda^2 || Rx||^2 } .
%
% Input variables
%  A -> Forward transfer matrix.
%  R -> Regularization matrix.
%  y -> observation vector
%  lambda -> Regularization coefficient, which determines the contribution
%            of the || Rx||^2 term to the solution.
%
% Output variables
% x_reg -> Regularized solution
% errNorm -> || A x_reg - y ||^2
% constNorm -> || Rx||^2
%
% Author: Onder Nazim Onak
%
    if (min(lambda) < 0)
      error('lambda must be > 0')
    end
    
    if(size(A,1) ~= size(y,1))
      error('Dimension mismatch for A matrix and y vector')
    end
    
    if(size(A,2) ~= size(R,2))
      error('Dimension mismatch for A and R matrices')
    end
    
    [m n] = size(A);

    if(m>n)
        k = 0;
    else
        k = n-m ;
    end    
    
    [U,V,W,C,M] = gsvd(A,R);
    X = inv(W');
    
    x_reg = zeros(size(A,2),length(lambda));
    errNorm = zeros(length(lambda),1);
    constNorm = zeros(length(lambda),1);
    
    for i=1:length(lambda)
        [ x_reg(:,i)] = solve(U,V,W,C,M,X,y,lambda(i),k);
        errNorm(i,1) = norm(y-A*x_reg(:,i));
        constNorm(i,1) = norm(R*x_reg(:,i));
    end
end

%
% solve l2l2 regularization problem
%
function [x_sol] = solve(U,V,W,C,M,X,bF,lambda,k)

    lambda2 = lambda^2;
    
    sigmaA = sqrt(diag(C'*C));
    sigmaA = sigmaA(k+1:end);
    sigmaL = sqrt(diag(M'*M));
    sigmaL = sigmaL(k+1:end);

    beta = U'*bF;
    zeta = sigmaA.*beta;
    xi = zeta./(sigmaA.^2 + lambda2*sigmaL.^2);
    x_sol = X(:,k+1:end)*xi;
end
