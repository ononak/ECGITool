
function [x_est, ern,  cn] = LpLqReg2(A,L,b,p,q,lambda)
% LPLQREG2
% Solves the regularization problem
%      (1/p)|| Ax - y||p + (lambda^2/q)*||Lx||q
%
% Input: 
% A (m x n) forward matrix, 
% L (n x n) regularization matrix 
% p and q are the norms
% lambda is the regularization parameter.
%
% Output : 
% x_est Estimated x (nx1) vector
% ern  || Ax - y||p
% cn   ||Lx||q
%
% Reference: A GENERALIZED KRYLOV SUBSPACE METHOD FOR Lp-Lq MINIMIZATION
% A. LANZA, S. MORIGI, L. REICHEL, AND F. SGALLARI
% DOI:10.1137/140967982
%
%
% Author: Onder Nazim Onak ononak@gmail.com
%

lambda2 = lambda^2;

    if((min([p q]) <= 0) || (max([p q]) >2))
        error('p and q must be 0<p,q<=2');
    end
    
% init loop control params    
tol = 1e-3;
diff = 1e3;
MAX_ITER = 1000;
k = 0;

    % initial estimation
    x_est = (A'*A + lambda2*L'*L)\A'*b;
    % starting from initial estimation iterate till converge
  while (diff > tol) 
      
    v = abs(A*x_est - b);
    z = abs(L*x_est)+eps;

    wf = power(v,p-2);
    wr = power(z,q-2);
    WF = diag(wf);
    WR = diag(wr);

    lhs = A'*WF*A + lambda2*L'*WR*L;
    rhs = A'*WF*b;
    
    M2 = sqrt(WF);
    gsvd(A,L);
    x_est_new = lhs\rhs;
    
    difftmp = norm(x_est_new - x_est)/norm(x_est);
    
    if((k > MAX_ITER) || (abs(difftmp - diff) < tol))
        break;
    end
        
    k = k + 1;
    x_est = x_est_new;
       
    diff = difftmp; 
  end
  
  % error norm
  
  ern = norm(A*x_est - b,p);
  %constraint norm
  cn = norm(L*x_est,q);
end


