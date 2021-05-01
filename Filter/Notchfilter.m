function [fsignal] = Notchfilter(signal, fs,varargin)
% BWDFILTER Baseline Wander drift filter 
%  highpass filter (>0.5 Hz) designed to remove Baseline wander 
%
% Author: Onder Nazim Onak
%
% Input variables
% signal -> noisy signal 
% fstop -> stop band frequency
% fs -> sampling frequency
% ftype -> filter type: 
%  'ellip' for elliptic filter
%  'butter' for butterworth filter
%
% Optinal variables
% fstop -> stop band frequency (default 50 Hz)
% ftype -> Filter type 
%          'butter' -> Butterworth filter (default)
%          'ellip' -> Elliptic filter
%
%
% Output variables
% fsignal -> filtered signal

%%
switch nargin
    case 2
        fstop = 50;
        ftype = 'butterworth';        
    case 3  
        fstop = varargin{1};
        ftype = 'butterworth';
    case 4
        fstop = varargin{1};
        ftype = varargin{2};
    otherwise
        error('Unknown filter type');
end
    
%% Generate filter coefficients
num = [];
den = [];

fn = fs/2;
Ws = [(fstop - 1)  (fstop + 1)]/fs;                            % Normalised Stopband (Passband = fstop-1 Hz To fstop+1 Hz)

if(strcmp(ftype,'elliptic'))
    % elliptic filter
    Rp =  1;                                                    % Passband Ripple/Attenuation
    Rs = 50;                                                    % Stopband Ripple/Attenuation
    Wp = [(fstop - 4)  (fstop + 4)]/fs;                         % Normalised Passband (Passband = fstop-4 Hz To fstop+4 Hz)
    [n,Wp] = ellipord(Wp, Ws, Rp, Rs);                          % Calculate Elliptic Filter Optimum Order
    [z,p,k] = ellip(n, Rp, Rs, Wp,'stop');                      % Elliptic Filter
    [num,den] = zp2tf(z,p,k); 
    %freqz(num,den)
elseif(strcmp(ftype,'butterworth'))
     %Butterworth filter
    [num,den] = butter(3,Ws,'stop');
    %freqz(num,den)
else
    error('Unknown filter type')
end


%% Apply filter reverse and forward direction to obtain zero phase filtered signal

[fsignal] = filtfilt(num, den, signal);
end

