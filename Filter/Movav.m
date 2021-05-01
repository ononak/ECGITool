function [fsignal] = Movav(signal,window)
% MOVAV  one dimensional moving average filter
%
% Author: Onder Nazim Onak
%
% Input variables
% signal -> noisy signal of size
% window -> moving window size
%
% Output variables
% fsignal -> filtered signal

%% filter parameters
num = ones(1,window)./window;
den = [1];

%% Apply filter reverse and forward direction to obtain zero phase filtered signal
[fsignal] = filtfilt(num, den, signal);
end

