
function [time_delay] = ActivationTimeDelay(pot1, pot2)
% ACTIVATIONTIMEDELAY   Find the activation time differences between two
%                       leads based on cross correlation function. 
%
% Usage: 
%        time_delay = ActivationTimeDelay(pot1, pot2)
%
% Inputs:
%   pot1     First lead potential vector
%   pot2     Second lead potential vector
%
% Output:
%   time_delay  Activation delay between leads
%
% References       A method for determining high-resolution activation 
%                  time delays in unipolar cardiac mapping. Shors et al.
%                  DOI: ﻿10.1109/10.544343
%
%
% Author: Onder Nazim Onak ononak@gmail.com
%
    n = length(pot1);
    t = 1:1:n;

     
     dpot1 = DiffByRegularization(pot1);
     dpot2 = DiffByRegularization(pot2);
     
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

