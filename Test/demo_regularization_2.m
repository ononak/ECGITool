clear
clc
addpath('../GeometryTools')
addpath('../Regularization')
addpath('../ActivationRecoveryTime')
addpath('../Util')

% load utah data
load('utah.mat')



I = eye(size(ftm,2));


[ reg_param ] = GenerateRegParams(ftm , 50);
reg_param(reg_param < 1e-4)=[];

for time = 1:87

tic
[x_1,rho1,eta1] = Tikhonov(ftm,I,bspm(:,time),reg_param);
toc


[ ~, ~, regparamIndex1] = LCurveCorner(rho1, eta1, reg_param,false);

selected_reg_param(time) = reg_param(regparamIndex1);

x_sol(:,time) = x_1(:,regparamIndex1);
end


plot(selected_reg_param);




