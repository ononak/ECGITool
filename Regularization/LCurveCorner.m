function [ corner regparam regparamIndex] = LCurveCorner(errNorm, constrNorm, lambda, pltotion)
% LCURVECORNER Finds the corner of the L curve
%   
% Input variables
% errNorm -> list of || A x_reg - y ||^p 
% constNorm -> list of || Rx||^q
% lambda -> list of regularization coefficient
%
% Optinal variables
% pltotion -> boolean ploting option. If true plot ecg sinal and detected QRS
% points
%
%
% Output variables
% corner -> curvature
% reg_param -> Regularization parameter
% regparamIndex -> index of the reg_param in the lamda list
%
% Reference: The Use of the L-Curve in the Regularization of Discrete Ill-Posed Problem
% Per Christian Hansen and Dianne Prost O’Leary
% SIAM J. Sci. Comput., 14(6), 1487–1503
%
% Author: Onder Nazim Onak ononak@gmail.com

%% parameters
	if ~exist('pltotion')
		pltotion = false;
    end

    initial_lambda = lambda; 
    
    logErrNorm = log10(errNorm);
    logConstrNorm = log10(constrNorm);
    logErrNorm = SmoothingFilter( logErrNorm',3,11)';
    logConstrNorm = SmoothingFilter( logConstrNorm',3,11)';    
    
   %remove the repetative data values
   [logErrNorm,IA,IC] = unique(logErrNorm);
   logConstrNorm = logConstrNorm(IA);
   lambda = lambda(IA);
   [logConstrNorm,IA,IC] = unique(logConstrNorm);
    logErrNorm = logErrNorm(IA);
   lambda = lambda(IA);  
   
   N = 4;
   knots = zeros(N,length(lambda)-1);
   for i=1:length(lambda)-1     
     knots(:,i) = linspace(lambda(i),lambda(i+1),N);         
   end
   
   knots = unique(knots);

   %% fit picewise polynomials for residual and solution norms
   ppE = spline(lambda,logErrNorm);
   ppC = spline(lambda,logConstrNorm);
   
   ppEval = ppval(ppE,knots);
   ppCval = ppval(ppC,knots);

   %% compute derivative of the polynomials
   dppE = PolyDerivative(ppE);
   dppC = PolyDerivative(ppC);
   
   ddppE = PolyDerivative(dppE);
   ddppC = PolyDerivative(dppC);
   
   %% compute the curvature
   drho = ppval(dppE,knots);
   ddrho = ppval(ddppE,knots);
   deta = ppval(dppC,knots);
   ddeta = ppval(ddppC,knots);
   
   kappa = (drho.*ddeta - ddrho.*deta)./power((drho.^2 + deta.^2),3/2);
   
   % determine max curvature point and corresponding lamda
   [~, index] = max(kappa);
   corner = kappa(index(1));
   params = find(lambda<knots(index(1)));
   
   if(pltotion)
    plot(ppEval,ppCval);
   end
   
   if(isempty(params))
       warning("L-curve corner could not be determined");
       corner = [];
       regparam = [];
       regparamIndex = [];
       return;
   end
   
   regparam = lambda(params(1));
   regparamIndex = params(1);
    
   if(pltotion)
        line([logErrNorm(regparamIndex);logErrNorm(regparamIndex)], [logConstrNorm(1);logConstrNorm(end)],'linestyle','--')
        line([logErrNorm(1);logErrNorm(end)], [logConstrNorm(regparamIndex);logConstrNorm(regparamIndex)],'linestyle','--')
        line(logErrNorm, logConstrNorm, 'Marker', 'x', 'Color', 'r')

        txt = strcat('\lambda = ',sprintf('%.2f',regparam));
        text(logErrNorm(regparamIndex),logConstrNorm(end),txt)
   end
   
    pi = find(initial_lambda==regparam);
    regparamIndex = pi(1);
end


