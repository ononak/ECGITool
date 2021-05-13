
function [xk] = TikhLp(A,b, x0,p,lambda)
% TIKHLP Solves the regularization problem
%      || Ax - b||2 + lambda*||x||p
%
% Input: A (m x n) forward matrix, p is the norm and lambda iste regularization parameter
% Output : Estimated x (nx1) vector
%
% Reference: A projection-based algorithm for L2-Lp Tikhonov
% regularization, Leonardo S. Borges Fermín S.V. Bazán Luciano Bedin
% DOI: 10.1002/mma.5110
%
% Author: Onder Nazim Onak

tol = 1e-3;
res = 1e3;
k = 0;
MaxIter = 1500;
c = 1.1;
tau = 1.1;


% [m n] = size(A);
 [U, s, V] = svd(A);

xk = x0;

while (k < MaxIter)
    xtk = xestimate(A,U,s,V, c, b, lambda, p, xk);
    fxt = objective(b,A,xtk,lambda,p);

    fx = objective(b,A,xk,lambda,p);

    cc = c;
    while(fxt > fx)
       cc = tau*cc; 
       xtk = xestimate(A,U,s,V, cc, b, lambda, p, xk);
       fxt = objective(b,A,xk,lambda,p);
       if(c >1e100)
           break;
       end
    end

    dif(k+1) = norm(xk-xtk);
    den = norm(xtk);
    if(dif(k+1) == 0)
        break;
    end
    xk = xtk;

    if(k > 1)
    if( dif(k)/den < tol)
        break;
    end
    end

k = k+1;
end

end

function [fx] = objective(b, A, x, lambda, p)

     fx = norm(b - A*x)^2 + (lambda^2)*power(norm(x,p),p);       
end

function [xhat] = xestimate(A,U,S,V, c, b, lambda, p, xprevious)


 n = length(xprevious);
Mc = A'*A+c*eye(n);
 S=S'*S;
 Sc = eye(n)/c;
   for i = 1:n
       Sc(i,i) = 1/(S(i,i)+c);
   end
 %Sctmp = 1./diag(s'*s + c*eye(192));
 %Sc = diag(Sctmp);
 %Mci = V*Sc*V';
 Mci = inv(Mc);
 Dxdiag = power(abs(xprevious),p/2 -1);
 
 Dxdiag(Dxdiag == 0) = 1;
 Dx = diag(Dxdiag);
 DDx = Dx'*Dx;
 
 xhat = Mci*(A'*b + (c*eye(n) - 0.5*(p*(lambda^2))*DDx)*xprevious);
end
 

