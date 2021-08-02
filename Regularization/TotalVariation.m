function [ x_reg, errNorm, constNorm ] = TotalVariation( A, R, y, lambda)
% TOTALVARIATION
%  Solves the Total Variation regularization problem
%  min { || A x - y ||^2 + lambda^2 || Rx||^1 } .
%
% Input variables
%  A -> Forward transfer matrix.
%  R -> Regularization matrix. (First or second order derivative operator)
%  y -> observation vector
%  lambda -> Regularization coefficient, which determines the contribution
%            of the || Rx||^1 term to the solution.
%
% Output variables
% x_reg -> Regularized solution
% errNorm -> || A x_reg - y ||^2
% constNorm -> || Rx||^1
%
% Author: Onder Nazim Onak

   if (min(lambda) < 0)
      error('lambda must be > 0')
    end
    
    if(size(A,1) ~= size(y,1))
      error('Dimension mismatch for A matrix and y vector')
    end
    
    if(size(A,2) ~= size(R,2))
      error('Dimension mismatch for A and R matrices')
    end
    
    x_reg = zeros(size(A,2),length(lambda));
    errNorm = zeros(length(lambda),1);
    constNorm = zeros(length(lambda),1);
    
    for (i=1:length(lambda))
        [ x_reg(:,i), errNorm(i,1), constNorm(i,1)] = LpLqReg(A,R,y,2,1,lambda(i));
    end
end

