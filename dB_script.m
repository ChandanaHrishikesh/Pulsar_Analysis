clear all

% run_driver_simpuls_Tobs_All_new1
% driver to run the analisys:

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

addpath(genpath('/media/cta/EXTERNAL_USB/Chandana_Work/Snag/Snag'))

CRThr = 2.5;
day2sec = 86400;
Tp = 1; % time window in days
dt = 1/64; % sampling time in days^-1
%SKY  position and spot on the sky (ROI) centered on the ra and dec:

ra = '06:33:44.2159045'; % hour
dec = [06 32 34.94512]; % deg
hour2deg(ra)
dec(1) + dec(2)/60 + dec(3)/3600

foldername = sprintf("dB_%s%s%s%s_%02d%02d", ra(1), ra(2), ra(4), ra(5), dec(1), dec(2))
mkdir(foldername)

mp_ra = 0;
mp_dec = 0;
ROI = 0.1; % spot on the sky
n_pos = 2; % number of (linear) sky positions

% Time interval of the analisys:
T_start = 54682.655439529444266;
T_stop = 54700.0;%59624.61388888889;
T_obs = (T_stop - T_start) * day2sec; % observing time in seconds

% Minimum energy to select the photons:
E_min = 550;

% Source parameters:
% reference epoch of the frequency source:
PEPOCH = 54945;

delta_f0dot = 1/(2*T_obs*Tp*86400);

%known source parameters:
% f0 = 3.3625291587551729515;
% spin_down = -8.999148e-13;

delta_lambda0 = (delta_f0dot*2*dt)/2;
% frband = [3.362 3.363]; % frequency band to be analised
% f_min = frband(1);% minimum frequency can be changed
% f_max = frband(2);% max frequ of the spectra
% Nspin = 1;
% spin_down_interval = [spin_down-delta_f0dot*Nspin spin_down+Nspin*delta_f0dot];
% if spin_down_interval(1)<0
%    lambda_min=(spin_down_interval(1))/(2*f_min);
% else
%    lambda_min=(spin_down_interval(1))/(2*f_max);
% end
% 
% if spin_down_interval(2)<0
%    lambda_max=(spin_down_interval(2))/(2*f_max);
% else
%    lambda_max=(spin_down_interval(2))/(2*f_min);
% end
lambda_min = -7.0e-13;
lambda_max = -7.0e-14;
lambda_par = [lambda_min lambda_max];
lambda_points = round((lambda_max-lambda_min)/delta_lambda0)

%% Chosing lambda_num

lambda_num = 7;

param.input.Tp = Tp;
param.input.ra = ra;
param.input.dec = dec;
param.input.ROI = ROI;
param.input.Tstart = T_start;
param.input.Tstop = T_stop;
param.input.E_min = E_min;
param.input.PEPOCH = PEPOCH;
param.input.n_pos = n_pos;
param.input.dt = dt;
param.input.lambda = lambda_par;
param.input.CRTH = CRThr;
param.input.lambda_num = lambda_num;
% PARAMETRI in INPUT (Input parameters)
%
% ra ascensione retta della sorgente in ore(ra right ascension of the source in hours)
% dec declinazione della sorgente in ore ma deve esser data come un vettore(dec declination of source in hours but must be given as a vector)
%
% ROI raggio,intorno in gradi usato per (ra,dec) (ROI radius, around in degrees used for (ra, dec))
%
% T_start tempo di inizio dei dati relativi alla pulsar (T_start start time of pulsar data)
% T_stop  tempo di fine dei dati relativi alla pulsar(T_stop end time of pulsar data)
%
% E_min energia minima (E_min minimum energy)

%% Counting photons

tic

('SELEZIONA I FOTONI DELLA SORGENTE')
%Using source_photon function
[t_fotoni_mjd, ra_deg, dec_deg] = source_photon(ra, dec, mp_ra, mp_dec, PEPOCH, ROI, T_start, T_stop, E_min); % t dei fotoni che vengono da una fissata direzione (t of photons coming from a fixed direction)
param.t_fotoni_mjd = t_fotoni_mjd;
param.ra_deg = ra_deg;
param.dec_deg = dec_deg;

save('dB_current_RUN.mat', 'foldername', 'param', '-v7.3');

toc
