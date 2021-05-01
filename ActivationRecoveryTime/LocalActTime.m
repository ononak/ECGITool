function [ lat] = LocalActTime(signal)
%LOCALACTTIME Computes the local activation time
%
% Author: Onder Nazim Onak ononak@gmail.com

filteredSignal = zeros(size(signal));
[nlead ntime] = size(signal);
lat = zeros(nlead,1);
 
 t = 1:1:size(signal,2);
 
     for i=1:nlead       
      [dif difHR fsignal] = DiffByRegularization(signal(i,:));
       h = ntime/(length(difHR(1,:)) - 1); 
      [~, ind] = min(difHR);
      [lat(i,1)] = ind*h;
     end

end


