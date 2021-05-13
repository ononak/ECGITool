% Author: Onder Nazim Onak
addpath('../Filter');
load('../Test/noisy_bspm.mat')
%% moving average filter
t = 0:0.11:20;
x = sin(t);  
n = rand(1,length(t));
x = x+n;
    
i= 2;
x = ts.potvals(i,1:3000);

ws = 15;
fsignal = Movav(x,ws);

figure(1);
plot([x' fsignal'])
title(['moving average filter' ' window size = ' num2str(ws)])
legend({'noisy' 'filtered'})

%% High pass filter for base line wander drift
i =2;
filterType = 'butterworth';
[fsignal] = Bwdfilter(ts.potvals(i,:),2048,filterType);
figure(2);
plot([ts.potvals(i,:)' fsignal'])

title(['High pass filter for base line wander drift' ' Filter type = ' filterType])
legend({'noisy' 'filtered'})

%% band stop filter

fs =2048

i =1;
filterType = 'butterworth';
noisy = ts.potvals(i,1:1000);
fstop = 50;
[fsignal] = Notchfilter(noisy, fs,fstop);

figure(3);
plot([noisy' fsignal'])

title([num2str(fstop) ' Hz Notch filter' ' Filter type = ' filterType])
legend({'noisy' 'filtered'})
