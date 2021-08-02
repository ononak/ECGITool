clear
clc

load('noisy_data_2.mat') 

% Plot the input signal.
figure(1)
clf
bookfonts
plot(t,m_true,'k');
xlabel('Time (s)');
ylabel('Seismic Reflector Amplitude');
title('True model')
disp('Displaying the true model (fig. 1)')


% Plot the noisy signal.
figure(2)
clf
bookfonts
plot(t,dn_measurement,'k');
xlabel('Time (s)');
ylabel('Seismic Amplitude');
title('Noisy measurements')
disp('Displaying the noisy data (fig. 2)')


n = length(dn_measurement);
L=eye(n);

%generate Regularization parameters
reg_param =10.^(-4:.05:4)';

% denoise signal with Tikhonov regularization

[x_1,rho1,eta1] = Tikhonov(G,L,dn_measurement,reg_param);

[ ~, ~, regparamIndex1] = LCurveCorner(rho1, eta1, reg_param,true);

figure(3)
clf
bookfonts
plot(t,x_1(:,regparamIndex1),'k');
xlabel('Time (s)');
ylabel('d(t)');
title('Estimated Seismic Reflector Amplitude by Tikhonov Method')

reg_param =10.^(-3:.125:-1.5);
[x_2,rho2,eta2] = TotalVariation(G,L,dn_measurement,reg_param);

[ ~, ~, regparamIndex2] = LCurveCorner(rho2, eta2, reg_param,true);

figure(4)
clf
bookfonts
plot(t,x_2(:,regparamIndex2),'k');
xlabel('Time (s)');
ylabel('d(t)');
title('Estimated Seismic Reflector Amplitude by L1 Regularization')





