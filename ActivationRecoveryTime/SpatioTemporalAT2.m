
function [at] = SpatioTemporalAT2(potvals, geom)
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


[nlead ntime] = size(potvals);

A = eye(ntime);
% initialize activation time vector
[dif difHR fsignal] = DiffByRegularization(potvals);


wdif = dif;
for i=1:ntime   
    [GF, GVal] = SurfaceGradientField(geom, potvals(:,i));
    wdif(:,i) = GVal.*wdif(:,i);
end


for i=1:nlead
    [~,at(i,1)] = min(wdif(i,:));
end


end

