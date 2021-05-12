function [earthModel] = DdModelingMars(train_ep, train_bsp)
% DdModelingMars Generate model function for EGM variables in terms of body
% surface leads based on MARS (Multivariate Adaptive Regression Splines).
%
% INPUT: 
%
% train_ep  -> (Ns x Ne) training data for EGM data, 
%  Ns is the number of samples, Ne is the number of EGM lead

% train_bsp  -> (Ns x Nb) training data for BSPM data, 
%  Ns is the number of samples, Nb is the number of body surface lead
%
%
% OUTPUT : 
%
% earthModel -> MARS model for each EGM
%
%
% Reference: A Novel Data-Adaptive Regression Framework Based on Multivariate Adaptive Regression
% Splines for Electrocardiographic Imaging
% Onder Nazim Onak, Taha Erenler, Yeşim Serinağaoğlu Doğrusöz
% IEEE TBME 
% Author: Onder Nazim Onak

if(~exist('earth'))  
    msg = 'Could not found Earth: Multivariate Adaptive Regression Splines (MARS) package.Add to path if it is available.Otherwise download and install Standalone C version for MATLAB from http://www.milbo.users.sonic.net/earth/';
    error(msg);
end


[~, nof_egm] = size(train_ep);

    earthModel = [earth()];

    disp('Modeling process started');
    
    for i =1:nof_egm
        earthModel(i) = earth();
        earthModel(i).nMaxDegree = 1; 
        earthModel(i).nMaxTerms = 192;
        %earthModel(i).Thresh = 0.000001;
        earthModel(i).Trace = 0;
        earthModel(i) = earthModel(i).train(train_bsp,train_ep(:,i));
        earthModel(i).bx = [];
        earthModel(i).Residuals =[];      
    end
    
    disp('Modeling process completed');
end

