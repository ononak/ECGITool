


%% Generate EGM model based on MARS
clc
clear
addpath('../DataDriven');
load('trainingData.mat')
load('utah.mat') % test data

%% generate model for EGM 10
earthModel_10 = DdModelingMars(training.Ep(:,10), training.Bsp);
% predict EGM 20 corresponding to from new Bspm
% BSPM must be (TxN) matrix where T is the number of time samples and N is the number of BSP lead
egm_10 = earthModel_10.predict(bspm');   

 subplot(2,1,1);
 plot([ep(10,:)' egm_10])
 title('EGM 10')
 legend('True' ,'Estimated');
 xlabel('time');
 ylabel('mV')
 
 
%% generate model for EGM 50
earthModel_50 = DdModelingMars(training.Ep(:,50), training.Bsp);
% predict EGM 50 corresponding to from new Bspm
egm_50 = earthModel_50.predict(bspm');  

 subplot(2,1,2);
 plot([ep(50,:)' egm_50])
  title('EGM 50')
 legend('True', 'Estimated');
  xlabel('time');
 ylabel('mV')
