
function [sparsity] = GiniMeasure(signal)
% GINIMEASURE Measures the sparsity of the given signal
%
% Input: 
% signal: one dimensional signal
% Output: 
% sparsity: value between 0 and 1. Sparsity increase if return value gets
% closer to 1.
%
% Reference:  Comparing Measures of Sparsity, Niall Hurley, Scott Rickard,
% IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 55, NO. 10, OCTOBER 2009
%
%
% Author: Onder Nazim Onak ononak@gmail.com
%

[M N] = size(signal);
if(M > N)
    signal = signal';
    N = M;
end
val_sorted = sort(abs(signal));
den = sum(val_sorted);
num = val_sorted.*((N - [1:1:N] + 0.5)/N);

sparsity = 1-2*sum(num/(den));
end