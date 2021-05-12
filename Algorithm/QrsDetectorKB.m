function  [qrs_blk] = QrsDetectorKB(ecg_signal,fs,pltotion)
%QRSDETECTORKB Knowledge based QRS detector
%
% Author: Onder Nazim Onak
% 
% INPUT 
% ecg_signal -> ecg signal
% fs -> sampling frequency
% pltotion -> ploting option - boolean - true if the user wants to plot result: 
%  'ellip' for elliptic filter
%  'butter' for butterworth filter
%
% Optinal variables
% pltotion -> boolean ploting option. If true plot ecg sinal and detected QRS
% points
%
%
% INPUT variables
% qrs_blk -> array of structure, contain location of Q,R, and S points for each detected
% beat.
%
%
% Reference: Fast QRS Detection with Optimized Knowledge-Based Method:
% Evaluation on 11 standart ECG Databases
% Muhamed Elgendi
% PLoS One. 2013;8(9):e73557. Published 2013 Sep 16. doi:10.1371/journal.pone.0073557

%% parameters
	if ~exist('pltotion')
		pltotion = false;
	end

W1 = 97; %ms
W2 = 611; %ms
beta = 0.8; % offset fraction

W1n = round(W1*1000/fs); % number of samples
W2n = round(W2*1000/fs);
W1n=35
W2n = 220

%% bandpass filter 8-20 Hz
num = [];
den = [];

fn = fs/2; 

f1 = 8;
f2 = 20;


[num den] = butter(3,f1/fn,'high');
ecg_h = filtfilt(num, den, ecg_signal);
%ecg_h = ecg_h/max(ecg_h);

[num den] = butter(3,f2/fn,'low');
ecg_f = filtfilt(num, den, ecg_h);


%% Algorithm
ecg_squared = ecg_f.^2;

num = ones(1,W1n)/W1n;
MAqrs = filtfilt(num, 1, ecg_squared);


num = ones(1,W2n)/W2n;
MAbeat = filtfilt(num, 1, ecg_squared);


z = mean(ecg_squared);
alpha = beta*z;
THR1 = MAbeat + alpha;
THR2 = W1n;
%plot([ecg_squared MAqrs MAbeat THR1])

BloksofInterest = zeros(1,length(MAqrs));

    for n =1: length(MAqrs)
        if MAqrs(n) > THR1
            BloksofInterest(n) = 1;
        end
    end
    
 % if BloksofInterest value change from 0 to 1 it is block starting point  
 Block_start = find(diff(BloksofInterest) == 1);  
 % if BloksofInterest value change from 1 to 0 it is block en point  
 Block_end = find(diff(BloksofInterest) == -1);
 
 if(Block_start(end) > Block_end(end))    
     Block_start = Block_start(1:end-1);
 end
 
 if(Block_start(1) > Block_end(1))    
     Block_end = Block_start(2:end);
 end

 block_sf_points = sort([Block_start Block_end]);
  
 qrs_index = 0;
 qrs_blk = {};
 for j =1:2:length(block_sf_points)

     block_len = block_sf_points(j +1) - block_sf_points(j);
     
     if (block_len >= THR2)
         qrs_index = qrs_index +1;
         qrs_blk{qrs_index}.Q = block_sf_points(j); %Q point
         qrs_blk{qrs_index}.S = block_sf_points(j + 1); % S point
         
         [~,ind] = max(ecg_squared(block_sf_points(j):block_sf_points(j+1)));
         qrs_blk{qrs_index}.R = block_sf_points(j) + ind - 1; % R point     
    end

 end
 
%% plot results
    if(pltotion)
        plot(ecg_signal);
        title('ECG data and Q,R,S points for detected beats')
        hold

         for k =1:qrs_index
             line([qrs_blk{k}.Q;qrs_blk{k}.Q], [min(ecg_signal);max(ecg_signal)],'linestyle','--','Color', 'r') 
             line([qrs_blk{k}.S;qrs_blk{k}.S], [min(ecg_signal);max(ecg_signal)],'linestyle','--','Color', 'k')
             line(qrs_blk{k}.R, ecg_signal(qrs_blk{k}.R),'Marker','*','Color', 'm')

         end
    end
 
end

