function [fsignal] = Notchfilter(signal, fs, fstop)
% Notchfilter Notch filter
%   
%
% Author: Onder Nazim Onak
%
% Input variables
% signal -> noisy signal 
% fs -> sampling frequency
% fstop -> stop band frequency
%
%
% Output variables
% fsignal -> filtered signal

 
%% Generate filter coefficients
fn = fs/2; % Nyquist frequency
Ws = fstop/fn; % Normalised Stopband
bw = (fstop + 1)/fn; 
   
[num,den]=iirnotch(Ws,bw);
%fvtool(num,den)
%% Apply filter reverse and forward direction to obtain zero phase filtered signal

[fsignal] = filtfilt(num, den, signal);
end

