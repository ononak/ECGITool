function [ecg_blk] = EcgDelineator(ecg_signal,fs,pltotion)
%ECGDELINEATOR Knowledge based QRS detector
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
% Reference 1: Fast QRS Detection with Optimized Knowledge-Based Method:
% Evaluation on 11 standart ECG Databases
% Muhamed Elgendi
% PLoS One. 2013;8(9):e73557. Published 2013 Sep 16. doi:10.1371/journal.pone.0073557
%
% Reference 2: A Proof-of-Concept Study: Simple and Effective Detection of P and T Waves 
% in Arrhythmic ECG Signals
% Mohamed Elgendi, Marianna Meo and Derek Abbott 
% Bioengineering 2016, 3(4), 26; https://doi.org/10.3390/bioengineering3040026

%% parameters
	if ~exist('pltotion')
		pltotion = false;
    end
    
interval_before_R_peak_ms = 83; %msec
interval_after_R_peak_ms = 166; %msec

interval_before_R_peak_ns = round((fs * interval_before_R_peak_ms)/1000);
interval_after_R_peak_ns = round((fs * interval_after_R_peak_ms)/1000);

W1 = 55; %msec half of the P wave duration
W2 = 110; %msec P ave duration

W1n = round((fs * W1)/1000); % number of sample
W2n = round((fs * W2)/1000); % number of sample

PminR = 55; %msec
PmaxR = 470; %msec

PminRn = round((fs * PminR)/1000); % number of sample
PmaxRn = round((fs * PmaxR)/1000); % number of sample

RTmin = 110; %msec
RTmax = 860; % msec

RTminn = round((fs * RTmin)/1000); % number of sample
RTmaxn = round((fs * RTmax)/1000); % number of sample

%remove baseline
[f_ecg_signal] = Bwdfilter(ecg_signal,fs);
    
 %% First step: Detect Q R and S points and and remove them from signal to make the P and T wave dominant   
    [qrs_blk] = QrsDetectorKB(f_ecg_signal,fs,false);


    %% bandpass filter 0.5-10 Hz
fn = fs/2; 

f1 = 0.5;
f2 = 10;


[num, den] = butter(3,f1/fn,'high');
ecg_h = filtfilt(num, den, f_ecg_signal);
%ecg_h = ecg_h/max(ecg_h);

[num, den] = butter(3,f2/fn,'low');
ecg_f = filtfilt(num, den, ecg_h);

% Remove QRS block from signal to make the P and T wave dominant   
number_of_qrs_blk = length(qrs_blk);

for i = 1:number_of_qrs_blk   

  sp = qrs_blk{i}.R - interval_before_R_peak_ns;
  ep = qrs_blk{i}.R + interval_after_R_peak_ns;

  if(sp <= 0)
      sp = 1;
  end
  if(ep > length(ecg_f))
      ep = length(ecg_f);
  end

  ecg_f(sp:ep) = 0;  
end


%% Algorithm

num = ones(1,W1n)/W1n;
MAPeak = filtfilt(num, 1, ecg_f);


num = ones(1,W2n)/W2n;
MAPwave = filtfilt(num, 1, ecg_f);
 
% determine block of interest
BloksofInterest = zeros(1,length(MAPeak));

    for n =1: length(MAPeak)
        if MAPeak(n) > MAPwave(n)
            BloksofInterest(n) = 1;
        end
    end

 % P wave search   
  P_block = zeros(1,length(BloksofInterest));

  for j = 1:length(qrs_blk)
     for i = 1:length(BloksofInterest)
          PRdist = i - qrs_blk{j}.R;
          if ((PRdist <= -PminRn) && (PRdist >= -PmaxRn))
              P_block(i) = BloksofInterest(i);
          end
      end
  end

   % if BloksofInterest value change from 0 to 1 it is block starting point  
 P_block_start = find(diff(P_block) == 1);  
 % if BloksofInterest value change from 1 to 0 it is block en point  
 P_block_end = find(diff(P_block) == -1);    
 
 if(P_block_start(end) > P_block_end(end))    
     P_block_start = P_block_start(1:end-1);
 end
 
 if(P_block_start(1) > P_block_end(1))    
     P_block_end = P_block_start(2:end);
 end
 
 % P wave search intervals   
 P_block_sf_points = sort([P_block_start P_block_end]);
 p_index = 0;
 
 for i=1:2:length(P_block_sf_points)    
   p_block_len = P_block_sf_points(i +1) - P_block_sf_points(i);  
   if(p_block_len >= 0.75*W1n)
        p_index = p_index + 1;
        p_blk_tmp{p_index}.start = P_block_sf_points(i);
        p_blk_tmp{p_index}.end = P_block_sf_points(i +1);   
        [~,maxi] = max(MAPeak(p_blk_tmp{p_index}.start:p_blk_tmp{p_index}.end));
        p_blk_tmp{p_index}.peak = p_blk_tmp{p_index}.start + maxi;
   end   
 end
 
 
 % If more than one R wave detected in one RR interval, keep only one of
 % them having the highest peak value
 p_index = 0;
 for i = 1:length(qrs_blk) - 1   
   p_index = p_index + 1;
   current_peak = 0;      

     for j = 1:length(p_blk_tmp)
         
        if((qrs_blk{i}.R < p_blk_tmp{j}.peak) && (qrs_blk{i + 1}.R > p_blk_tmp{j}.peak))
            
             new_peak = f_ecg_signal(p_blk_tmp{j}.peak);     
             
             if(new_peak >= current_peak)
                p_blk{p_index} = p_blk_tmp{j};
                current_peak = f_ecg_signal(p_blk{p_index}.peak);
             end
        elseif(qrs_blk{i + 1}.R < p_blk_tmp{j}.peak)
            break;
        end
        
     end
 end 
 
 
 % Twave search
  
 T_block = zeros(1,length(BloksofInterest));

 for j = 1:length(qrs_blk)    
     for i = 1:length(BloksofInterest)
          RTdist = qrs_blk{j}.R - i;
          if ((RTdist <= -RTminn) && (RTdist >= -RTmaxn))
              T_block(i) = BloksofInterest(i);
          end
     end
 end 

  % if BloksofInterest value change from 0 to 1 it is block starting point  
 T_block_start = find(diff(T_block) == 1);  
 % if BloksofInterest value change from 1 to 0 it is block en point  
 T_block_end = find(diff(T_block) == -1);    
 
 if(T_block_start(end) > T_block_end(end))    
     T_block_start = T_block_start(1:end-1);
 end
 
 if(T_block_start(1) > T_block_end(1))    
     T_block_end = T_block_start(2:end);
 end
 
 % T wave search intervals   
 T_block_sf_points = sort([T_block_start T_block_end]);
 
 t_index = 0;
 for i=1:2:length(T_block_sf_points)
   t_block_len = T_block_sf_points(i +1) - T_block_sf_points(i);  
   if(t_block_len >= 1.25*W1n)    
       t_index = t_index + 1;
       t_blk_tmp{t_index}.start = T_block_sf_points(i);
       t_blk_tmp{t_index}.end = T_block_sf_points(i +1);   
       [~,maxi] = max(MAPeak(t_blk_tmp{t_index}.start:t_blk_tmp{t_index}.end));
       t_blk_tmp{t_index}.peak = t_blk_tmp{t_index}.start + maxi;
   end
 end   
    
 % If more than one T wave detected in one RR interval, keep only one of
 % them having the highest peak value
 t_index = 0;
 for i = 1:length(qrs_blk) - 1   
   t_index = t_index + 1;
   current_peak = 0;      

     for j = 1:length(t_blk_tmp)
         
        if((qrs_blk{i}.R < t_blk_tmp{j}.peak) && (qrs_blk{i + 1}.R > t_blk_tmp{j}.peak))
            
             new_peak = f_ecg_signal(t_blk_tmp{j}.peak);     
             
             if(new_peak >= current_peak)
                t_blk{t_index} = t_blk_tmp{j};
                current_peak = f_ecg_signal(t_blk{t_index}.peak);
             end
        elseif(qrs_blk{i + 1}.R < t_blk_tmp{j}.peak)
            break;
        end
        
     end
 end
 
 
 ecg_blk.qrs_blk = qrs_blk;
 ecg_blk.p_blk = p_blk;
 ecg_blk.t_blk = t_blk;
 
 %% plot results
    if(pltotion)
        plot(ecg_signal);

        title('ECG data and P,Q,R,S,T points for detected beats')
        hold

         for k =1:number_of_qrs_blk
             %line([ecg_blk.qrs_blk{k}.Q;ecg_blk.qrs_blk{k}.Q], [min(f_ecg_signal);max(f_ecg_signal)],'linestyle','--','Color', 'r') 
             %line([ecg_blk.qrs_blk{k}.S;ecg_blk.qrs_blk{k}.S], [min(f_ecg_signal);max(f_ecg_signal)],'linestyle','--','Color', 'r')
             line(ecg_blk.qrs_blk{k}.Q, ecg_signal(ecg_blk.qrs_blk{k}.Q),'Marker','*','Color', 'm','MarkerSize',12)
             line(ecg_blk.qrs_blk{k}.R, ecg_signal(ecg_blk.qrs_blk{k}.R),'Marker','*','Color', 'm','MarkerSize',12)
             line(ecg_blk.qrs_blk{k}.S, ecg_signal(ecg_blk.qrs_blk{k}.S),'Marker','*','Color', 'm','MarkerSize',12)
         end
         
         for k =1:p_index
             %line([ecg_blk.p_blk{k}.start;ecg_blk.p_blk{k}.start], [min(f_ecg_signal);max(f_ecg_signal)],'linestyle','-','Color', 'k') 
             %line([ecg_blk.p_blk{k}.end;ecg_blk.p_blk{k}.end], [min(f_ecg_signal);max(f_ecg_signal)],'linestyle','-','Color', 'k')
             line(ecg_blk.p_blk{k}.peak, ecg_signal(ecg_blk.p_blk{k}.peak),'Marker','+','Color', 'r','MarkerSize',12)
         end
         
         for k =1:t_index
             %line([ecg_blk.t_blk{k}.start;ecg_blk.t_blk{k}.start], [min(f_ecg_signal);max(f_ecg_signal)],'linestyle','-','Color', 'c') 
             %line([ecg_blk.t_blk{k}.end;ecg_blk.t_blk{k}.end], [min(f_ecg_signal);max(f_ecg_signal)],'linestyle','-','Color', 'c')
             line(ecg_blk.t_blk{k}.peak, ecg_signal(ecg_blk.t_blk{k}.peak),'Marker','o','Color', 'r','MarkerSize',12)
         end         
         
    end
 
end
