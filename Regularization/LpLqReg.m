
function [x_hat, ern,  cn] = LpLqReg(A,R,b,p,q,lambda,x0)
% LPLQREG
% Solves the regularization problem
%      (1/p)|| A*x - y||p + (lambda^2/q)*||R*x||q
%
% Input: 
% A (m x n) forward matrix, 
% L (n x n) regularization matrix 
% p and q are the norms
% lambda is the regularization parameter.
%
% Output : 
% x_est Estimated x (nx1) vector
% ern  || A*x_est - y||p
% cn   ||R*x||q
%
% Reference: A GENERALIZED KRYLOV SUBSPACE METHOD FOR Lp-Lq MINIMIZATION
% A. LANZA, S. MORIGI, L. REICHEL, AND F. SGALLARI
% DOI:10.1137/140967982
%
%
% Author: Onder Nazim Onak ononak@gmail.com
%
    if(lambda < 0)
        error('lambda must be >0');
    end

    if((min([p q]) <= 0) || (max([p q]) >2))
        error('p and q must be 0<p,q<=2');
    end
  
% init loop control params    
tol = 1e-3;
diff = 1e10;
MAX_ITER = 200;

    % initial estimation
    if(nargin == 7)
        x_est = x0;
    else
     x_est = solve(A,R,b,lambda);
    end

  iter = 0;
    %% starting from initial estimation iterate till converge
    
    hold off
  while (diff > tol) 

  %% Compute weight matrices    
    v = abs(A*x_est - b);
    z = abs(R*x_est)+eps;

    wf = v.^((p-2)/2);
    wr = z.^((q-2)/2);

    WF = diag(wf);
    WR = diag(wr);

    AF = WF*A;
    LR = WR*R;
    bF = WF*b;

    %% solve weighted regularization
    x_est_new = solve(AF,LR,bF,lambda);

    %% check convergence
    difftmp = norm(x_est_new - x_est)/norm(x_est);

    if((iter > MAX_ITER) || (difftmp < tol))
        break;
    end
    
%     if( (difftmp > 1.2*diff))
%         break
%     end
    
    iter = iter + 1;
    x_est = x_est_new;    
    diff = difftmp;
    
  end

  x_hat = x_est;
  % error norm
  ern = norm(A*x_est - b,p); 
  %constraint norm
  cn = norm(R*x_est,q);
  
 end

%
% solve l2l2 regularization problem
%
function [x_sol] = solve(AF,LR,bF,lambda)
   
    [m, n] = size(AF);

    if(m>n)
        k = 0;
    else
        k = n-m ;
    end
    
    lambda2 = lambda^2;
    
    [U,V,W,C,M] = gsvd(AF,LR);
    X = inv(W');

    sigmaA = sqrt(diag(C'*C));
    sigmaA = sigmaA(k+1:end);
    sigmaL = sqrt(diag(M'*M));
    sigmaL = sigmaL(k+1:end);

    beta = U'*bF;
    zeta = sigmaA.*beta;
    xi = zeta./(sigmaA.^2 + lambda2*sigmaL.^2);
    x_sol = X(:,k+1:end)*xi;

end

