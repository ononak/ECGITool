% Author: Onder Nazim Onak
clc
clear
addpath('../Filter');
load('../Test/noisy_bspm.mat')

fs = 2048;

%% moving average filter  
i= 2;
x = ts.potvals(i,:);

ws = 15;
fsignal = Movav(x,ws);

figure(1);
subplot(2,1,1)
plot([x' fsignal'])
title(['moving average filter' ' window size = ' num2str(ws)])
xlabel('msec');
ylabel('mV');
legend({'noisy' 'filtered'})

% power spectrum
ynoisy = fft(x);
yfiltered = fft(fsignal);
n = length(x);          % number of samples
fstep = (0:n-1)*(fs/n);     % frequency range
power_noisy = abs(ynoisy).^2/n;    % power of the DF
power_filtered = abs(yfiltered).^2/n;    % power of the DF

subplot(2,1,2)
plot(fstep(100:500),[power_noisy(100:500)' power_filtered(100:500)'],'LineWidth',2)
title('Power spectrum')
xlabel('Frequency');
ylabel('Power');
legend({'noisy' 'filtered'})

%% High pass filter for base line wander drift
i =2;
filterType = 'butterworth';
[fsignal] = Bwdfilter(ts.potvals(i,:),fs,filterType);



figure(2);
subplot(2,1,1)
plot([ts.potvals(i,:)' fsignal'])
title(['High pass filter for base line wander drift' ' Filter type = ' filterType])
xlabel('msec');
ylabel('mV');
legend({'noisy' 'filtered'})

% power spectrum
ynoisy = fft(ts.potvals(i,:));
yfiltered = fft(fsignal);
n = length(ts.potvals(i,:));          % number of samples
fstep = (0:n-1)*(fs/n);     % frequency range
power_noisy = abs(ynoisy).^2/n;    % power of the DF
power_filtered = abs(yfiltered).^2/n;    % power of the DF

subplot(2,1,2)
plot(fstep(1:50),[power_noisy(1:50)' power_filtered(1:50)'],'LineWidth',2)
title('Power spectrum')
xlabel('Frequency');
ylabel('Power');
legend({'noisy' 'filtered'})



%% band stop filter

i =1;
filterType = 'butterworth';
noisy = ts.potvals(i,:);
fstop = 50;
[fsignal] = Notchfilter(noisy, fs,fstop);

figure(3);
subplot(2,1,1)
plot([noisy' fsignal'])

title([num2str(fstop) ' Hz Notch filter' ' Filter type = ' filterType])
xlabel('msec');
ylabel('mV');
legend({'noisy' 'filtered'})


% power spectrum
ynoisy = fft(noisy);
yfiltered = fft(fsignal);
n = length(ts.potvals(i,:));          % number of samples
fstep = (0:n-1)*(fs/n);     % frequency range
power_noisy = abs(ynoisy).^2/n;    % power of the DF
power_filtered = abs(yfiltered).^2/n;    % power of the DF

subplot(2,1,2)
plot(fstep(100:350),[power_noisy(100:350)' power_filtered(100:350)'],'LineWidth',2)
title('Power spectrum')
xlabel('Frequency');
ylabel('Power');
legend({'noisy' 'filtered'})
