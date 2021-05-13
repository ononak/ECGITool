function [fsignal] = Bwdfilter(signal,fs,varargin)
% BWDFILTER Baseline Wander drift filter 
%  highpass filter (>0.5 Hz) designed to remove Baseline wander 
%
% Author: Onder Nazim Onak
% uu
% Input variables
% signal -> noisy signal 
% fs -> sampling frequency
% ftype -> filter type: 
%  'ellip' for elliptic filter
%  'butter' for butterworth filter
%
% Optinal variables
% ftype -> Filter type 
%          'butter' -> Butterworth filter
%          'ellip' -> Elliptic filter
%
%
% Output variables
% fsignal -> filtered signal

%%
switch nargin
    case 2    
        ftype = 'butterworth';
    case 3
        ftype = varargin{:};
    otherwise
        error('Unknown filter type');
end
    
%% Generate filter coefficients
num = [];
den = [];

fn = fs/2;                                                  % Nyquist Frequency (Hz)

if(strcmp(ftype,'elliptic'))
% Eliptic filter params cut off freq 0.5, passband attenuation 3 dB
% stopband attenuation 60dB
 [num, den] = ellip(3,3,60,.5/fn,'high');
elseif(strcmp(ftype,'butterworth'))
 %Butterworth filter
 [num, den] = butter(2,.5/fn,'high');
else
    error('Unknown filter type')
end
  
%% Apply filter reverse and forward direction to obtain zero phase filtered signal

[fsignal] = filtfilt(num, den, signal);

end

