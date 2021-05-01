clear
clc
addpath('../GeometryTools')
addpath('../Regularization')
addpath('../ActivationRecoveryTime')
addpath('../Util')

% load utah data
load('utah.mat')

time = 15;
I = eye(size(ftm,2));
D = SurfaceGradient(epigeom);
L = SurfaceLaplacian(epigeom);

[ reg_param ] = GenerateRegParams(ftm , 50);
reg_param(reg_param < 1e-4)=[];

fprintf('Step 1: Solve inverse problem using %d different regularization parameters\n',length(reg_param));
fprintf('Step 2: Select the best solution based on L-curve criteria\n\n');

disp('Zero Order Tikhonov')
tic
[x_1,rho1,eta1] = Tikhonov(ftm,I,bspm(:,time),reg_param);
toc
disp('First Order Tikhonov')
tic
[x_2,rho2,eta2] = Tikhonov(ftm,D,bspm(:,time),reg_param);
toc
disp('Second Order Tikhonov')
tic
[x_3,rho3,eta3] = Tikhonov(ftm,L,bspm(:,time),reg_param);
toc
disp('Total Variation with Surface Gradient operator')
tic
[x_4,rho4,eta4] = TotalVariation(ftm,D,bspm(:,time),reg_param);
toc
disp('Total Variation with Surface Laplace operator')
tic
[x_5,rho5,eta5] = TotalVariation(ftm,L,bspm(:,time),reg_param);
toc

[ ~, ~, regparamIndex1] = LCurveCorner(rho1, eta1, reg_param);
[ ~, ~, regparamIndex2] = LCurveCorner(rho2, eta2, reg_param);
[ ~, ~, regparamIndex3] = LCurveCorner(rho3, eta3, reg_param);
[ ~, ~, regparamIndex4] = LCurveCorner(rho4, eta4, reg_param);
[ ~, ~, regparamIndex5] = LCurveCorner(rho5, eta5, reg_param);

figure;
plot([ ep(:,time) x_1(:,regparamIndex1) x_2(:,regparamIndex2) x_3(:,regparamIndex3) x_4(:,regparamIndex4) x_5(:,regparamIndex5)]);
legend('True', 'Zot','Fot','Sot','TV1','TV2');
title('Estimation results')
ylabel('mV')
xlabel('lead')


[at_true] = SpCoherentActTime(ep, epigeom);
figure;
patch('Vertices',epigeom.pts','Faces',epigeom.fac','FaceVertexCData',at_true,'FaceColor','interp');
title('True AT map')

% [at_zot] = SpCoherentActTime(x_1, epigeom);
% figure;
% patch('Vertices',epigeom.pts','Faces',epigeom.fac','FaceVertexCData',at_zot,'FaceColor','interp');
% title('AT map obtained by using Zot estimation results')
% 
% [at_tv1] = SpCoherentActTime(x_4, epigeom);
% figure;
% patch('Vertices',epigeom.pts','Faces',epigeom.fac','FaceVertexCData',at_tv1,'FaceColor','interp');
% title('AT map obtained by using TV1 estimation results')

