clear
clc

load('noisy_data_1.mat') 

% Plot the input signal.
figure(1)
clf
bookfonts
plot(t,m_true,'k');
xlabel('Time (s)');
ylabel('m(t)');
title('True model')
disp('Displaying the true model (fig. 1)')


% Plot the noisy signal.
figure(2)
clf
bookfonts
plot(t,noisy_observation,'k');
xlabel('Time (s)');
ylabel('d(t)');
title('Noisy model')
disp('Displaying the noisy data (fig. 2)')


n = length(noisy_observation);
G=eye(n);

%generate Regularization parameters
reg_param =10.^(-2:.25:5)';

% denoise signal with Tikhonov regularization

[x_1,rho1,eta1] = Tikhonov(G,D1,noisy_observation,reg_param);

[ ~, ~, regparamIndex1] = LCurveCorner(rho1, eta1, reg_param,false);

figure(3)
clf
bookfonts
plot(t,x_1(:,regparamIndex1),'k');
xlabel('Time (s)');
ylabel('d(t)');
title('Denoised signed with 1th-order Tikhonov')


[x_2,rho2,eta2] = Tikhonov(G,D2,noisy_observation,reg_param);

[ ~, ~, regparamIndex2] = LCurveCorner(rho2, eta2, reg_param,true);

figure(4)
clf
bookfonts
plot(t,x_2(:,regparamIndex2),'k');
xlabel('Time (s)');
ylabel('d(t)');
title('Denoised signed with 2nd-order Tikhonov')





