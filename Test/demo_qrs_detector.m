% Author: Onder Nazim Onak
%%
addpath('../Algorithm');
load('ecg2.mat')

disp(['ECG data sampling frequency ' num2str(ecg_s.fs)]);

[qrs_blk] = QrsDetectorKB(ecg_s.pval,ecg_s.fs,true);

disp(['Number of detected beat ' num2str(length(qrs_blk))]);