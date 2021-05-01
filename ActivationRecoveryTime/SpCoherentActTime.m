
function [at] = SpCoherentActTime(potvals, geom)
%
% SPCOHERENTACTTIME   Find the spatially coherent activation times 
%
% Usage: 
%        at = SpCoherentActTime(potvals, geom);
%
% Inputs:
%   potvals     Time varying potential values for leads (size: (m x n) m: number of leads, n: number of time instance)
%   geom        Heart surface geometry. Must include triangle information geom.fac (size: (3xnumOfTri); 
%               node indexing should start with 1, not 0)
%
% Output:
%   at          Activation time vector (size: m x 1):
%
% References:      Spatially Coherent Activation Maps for Electrocardiographic Imaging, ﻿
%                  Duchateau et. al DOI:10.1109/TBME.2016.2593003
%
% Author: Onder Nazim Onak ononak@gmail.com
%

% take the derivative of the potentials

[m, ~] = size(potvals);

% filter potvals and compute time derivatives
[dif difHR fsignal] = DiffByRegularization(potvals);

    neighbours = NeighbourList(geom.fac);

    % For each measurement point calculate the activation time delays between
    % its first order neighbours
    [mn, ~] = size(neighbours);
    neigSize = sum(neighbours(:,2))/2;
    D_incidence = spalloc(neigSize,mn, neigSize*2);
    FON_time_delay = zeros(neigSize,1); % time delay vector
    index = 1;
    for i =1:mn
        for j = 1:neighbours(i,2)
            if(i < neighbours(i,j+2))
            D_incidence(index,i) = 1;
            D_incidence(index,neighbours(i,j+2)) = -1;
            FON_time_delay(index,1) = TimeDelay(dif(i,:),dif(neighbours(i,j+2),:));
              index = index +1;
            end
        end
    end

    % Find the local activation times
   LAT = Lat(dif, difHR, fsignal);
    
    % Recast the activation times considering the computed time delays 
    I = eye(m,m);
    
    A = [I; D_incidence];
    C = [LAT;FON_time_delay];
    
    at = A\C;
       
end



function [time_delay] = TimeDelay(dpot1, dpot2)
%    Compute the time delay between two
%    signal based on cross correlation function. 
%
%  dpot1 and dpot2 are signal time derivatives

    n = length(dpot1);

    % compute cross correlation for electrograms pair   
    [acor,lag] = xcorr(dpot1,dpot2);
    
    % compute Hilbert trasform of the cross correlation between the
    % electrograms
    %Hacor = hilbert(acor);
    
    % Use interpolation to increase the resolution 
    F = griddedInterpolant(lag,acor);
    xq = linspace(lag(1),lag(end),n*100);
    vq = F(xq);
    % Find maximum point of the cross correlation function
    [~, ind] = max(vq);
    time_delay = xq(ind);
end

function [ lat] = Lat(dif, difHR, fsignal)
%LOCALACTTIME Computes the local activation time

[nlead, ntime] = size(fsignal);
lat = zeros(nlead,1);
h = ntime/(length(difHR(1,:)) - 1);
     for i=1:nlead        
      [~, ind] = min(difHR(i,:));
      [lat(i,1)] = ind*h;
     end

end

