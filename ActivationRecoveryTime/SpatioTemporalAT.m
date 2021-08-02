
function [at] = SpatioTemporalAT(potvals, L)
%Estimates activation times based on spatio-temporal regularization.
% Usage: 
%        TActivation = ActivationTime(potvals, geom);
%
% Inputs:
%   potvals     Time varying potential values for leads (size: (m x n) m: number of leads, n: number of time instance)
%   L        Heart surface surface Laplacian matrix
%
% Output:
%   at  Activation time vector (size: m x 1):
%
% Reference: Erem, B., Coll-Font, J., Orellana, R. M., St'Ovicek, P., &
% Brooks, D. H. (2014). Using transmural regularization and dynamic mod-
% eling for noninvasive cardiac potential imaging of endocardial pacing with
% imprecise thoracic geometry. IEEE Transactions on Medical Imaging, 33(3),
% 726-738. http://doi.org/10.1109/TMI.2013.2295220

[nlead ntime] = size(potvals);


A = eye(nlead);
% initialize activation time vector
[ lat] = LocalActTime(potvals);

[reg_param] = GenerateRegParams(A , 250, L);

[gat_tmp,rho,eta] = Tikhonov(A,L,lat,reg_param);
[ ~, ~, regparamIndex] = LCurveCorner(rho, eta, reg_param,true);

if(isempty(regparamIndex))
    regparamIndex = length(reg_param);
end

at = gat_tmp(:,regparamIndex);

end

