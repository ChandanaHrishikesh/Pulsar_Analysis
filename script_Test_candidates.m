% run_driver_simpuls_Tobs_All_new1
%driver to run the analisys:

% around a spot on the sky:
%
% around a frequency source OR over all the freqnency band (1/(2Tp):
% spin_par(1)=f0 OR spin_par(1)=0;
% a source OR blind over a frequency band: 

% % 
% % 
% % SORGENTE J0633+0632
% % 
% % [ampiez_cand_max freq_cand_max f0dot_cand_max lambda_fin_max RA_cand_max DEC_cand_max]=simpuls_blind('06:33:44.2159045',[06 32 34.94512],0,0,0.6,54682.655439529444266,55208.116494250447431,550,54945,[3.3625291587551729515 -8.999148066954602476e-13 -2.4123046863284875109e-23],6,10,7,1/64)
% % 
% % 
% % T_obs = 4.554619750376010e+07;
% % f0 = 3.362549545880845;
% % f1dot = -8.993681449405110e-13;

addpath(genpath('/data/virgo/Work_dantonio/Snag/'))
CRThr=10;
day2sec=86400;
Tp=15; % time window in days
dt=1/64; % sampling time in days^-1
%SKY  position and spot on the sky (ROI) centered on the ra and dec:

ra='06:33:44.2159045' % hour
dec=[06 32 34.94512]; % deg
hour2deg(ra)
dec(1)+dec(2)/60+dec(3)/3600

mp_ra=0
mp_dec=0
ROI=0.6 % spot on the sky
n_pos=0; % number of ...=6

% Time interval of the analisys:
T_start=54682.655439529444266;
T_stop=59624.61388888889;
T_obs=(T_stop-T_start)*day2sec %  observing time in seconds

% Minimum energy to select the photons:
E_min=550;

% Source parameters:
% reference epoch of the frequency source:
PEPOCH=54945;

delta_f0dot = 1/(2* T_obs*Tp*86400);

% KNOWN SOURCE FREQUENCY AND SPIN_DOWN: 

directed=0;
lambda_par=[0 0]; % not used
%known source parameters:  
f0=3.3625291587551729515 
spin_down=-8.999148e-13
spin_par(1)=f0;
spin_par(2)=spin_down;
spin_par(3)=0;
% "DIRECTED SEARCH":
directed=1;

% spin_down_interval=[-8.999148066954602476e-13-delta_f0dot*5 -8.999148066954602476e-13+delta_f0dot*5] 
% n_f0dot=round((spin_down_interval(2)-spin_down_interval(1))/delta_f0dot/2); % HALF  number of spin_down point

%% Evaluate number of lambda points = (lambda vero of a given source is: spin_down/(2f0) 
% f0 is the frequency at the PEPOCH of the source
spin_par=[0 0 0] % not used
delta_lambda=(delta_f0dot*2*dt)/2;
frband=[3 5]; % frequency band to be analised
f_min=frband(1);% minimum frequency can be changed
f_max=frband(2);% max frequ of the spectra
Nspin=5;
spin_down_interval=[spin_down-delta_f0dot*Nspin spin_down+Nspin*delta_f0dot];
if spin_down_interval(1)<0
    lambda_min=(spin_down_interval(1))/(2*f_min);
else
    lambda_min=(spin_down_interval(1))/(2*f_max);
end

if spin_down_interval(2)<0
    lambda_max=(spin_down_interval(2))/(2*f_max);
else
    lambda_max=(spin_down_interval(2))/(2*f_min);
end
lambda_par=[lambda_min lambda_max];
lambda_points=round((lambda_max-lambda_min)/delta_lambda)

%[ampiez_cand_max freq_cand_max f0dot_cand_max lambda_fin_max RA_cand_max DEC_cand_max]=simpuls_blind_spectrumTest(ra,dec,mp_ra,mp_dec,ROI,T_start,T_stop,E_min,PEPOCH,spin_par,n_pos,Tp,dt,lambda_par,allsky)
[ampiez_cand_max freq_cand_max f0dot_cand_max lambda_fin_max RA_cand_max DEC_cand_max]=simpuls_blind_Test_candidates(ra,dec,mp_ra,mp_dec,ROI,T_start,T_stop,E_min,PEPOCH,spin_par,n_pos,Tp,dt,lambda_par,directed,CRThr,spin_down,frband,T_obs,lambda_points)
  
%[ampiez_cand_max freq_cand_max f0dot_cand_max lambda_fin_max RA_cand_max DEC_cand_max]=simpuls_blind_spectrumTest(ra,dec,mp_ra,mp_dec,ROI,T_start,T_stop,E_min,PEPOCH,spin_par,n_pos,n_f0dot,Tp,dt)
%function [cl_grid x_grid ampiez_cand freq_cand f0dot_cand lambda_fin RA_cand DEC_cand]=simpuls_blind(ra,dec,mp_ra,mp_dec,ROI,T_start,T_stop,E_min,PEPOCH,spin_par,n_pos,n_f0dot,Tp,dt)

% SORGENTE J0633+0632:
% f0dot=-8.999148066954602476  10?^-13
% freq=3.3625291587551729515
%[ampiez_cand freq_cand f0dot_cand lambda_fin RA_cand DEC_cand]=simpuls_blind_spectrum('06:33:44.2159045',[06 32 34.94512],0,0,0.6,54682.655439529444266,59624.61388888889,550,54945,[0 -9.0e-13 0],6,10,7,1/64)
% 
% % SORGENTE J0633+0632
% POS: ('06:33:44.2159045',[06 32 34.94512],0,0,0.6,54682.655439529444266,59624.61388888889,550,54945,[3.3625291587551729515 -8.999148066954602476e-13 0],6,0,7,1/64,1)
% 
% [ampiez_cand freq_cand f0dot_cand lambda_fin RA_cand DEC_cand]=simpuls_blind_spectrumTest('06:33:44.2159045',[06 32 34.94512],0,0,0.6,54682.655439529444266,59624.61388888889,550,54945,[3.3625291587551729515 -8.999148066954602476e-13 0],6,0,7,1/64,1)
