
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

[m, n] = size(potvals);

[filteredSignal] = SmoothingFilter(potvals, 3);

dif = zeros(m,n);
for i=1:m
    dif(i,:) = gradient(filteredSignal(i,:));
end
    

wdif = dif;
for i=1:n   
    [GFV GFN] = SurfaceGradientField(geom, potvals(:,i));
    wdif(:,i) = GFN'.*dif(:,i);
end


for i=1:m
    [~,at(i,1)] = min(wdif(i,:));
end


end

