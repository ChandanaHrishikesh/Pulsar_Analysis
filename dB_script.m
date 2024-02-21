clear all

addpath(genpath('/data/virgo/Matlab/Snag/'))

currentfile = 'dB_current_RUN.mat'

CRThr = 4;
Tp = 7; % time window in days
dt = 1/64; % sampling time in days^-1

%S KY  position and spot on the sky (ROI) centered on the ra and dec:
ra = '06:33:44.2159045'; % hour
dec = [06 32 34.94512]; % deg
hour2deg(ra)
dec(1) + dec(2)/60 + dec(3)/3600

ROI = 0.6; % spot on the sky
n_pos = 3; % number of (linear) sky positions

% Time interval of the analisys:
T_start = 54682.655439529444266;
T_stop = 59624.61388888889;

% Minimum energy to select the photons:
E_min = 550;

% Source parameters:
% reference epoch of the frequency source:
PEPOCH = 54945;

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
param.input.CRTH = CRThr;
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

foldername = sprintf("dB_%s%s%s%s_%02d%02d", ra(1), ra(2), ra(4), ra(5), dec(1), dec(2))
mkdir(foldername);

save(currentfile, 'param', 'foldername', '-v7.3');

%% Estimating lambda_points

load(currentfile, 'param');

T_obs = (param.input.Tstop - param.input.Tstart) * 86400; % observing time in seconds

delta_f0dot = 1/(2 * T_obs * param.input.Tp*86400);

%known source parameters:
% f0 = 3.3625291587551729515;
 spin_down = -8.999148e-13;

delta_lambda0 = (delta_f0dot*2*param.input.dt)/2;
frband = [3.362 3.363]; % frequency band to be analised
f_min = frband(1);% minimum frequency can be changed
f_max = frband(2);% max frequ of the spectra
Nspin = 1;
spin_down_interval = [spin_down-delta_f0dot*Nspin spin_down+Nspin*delta_f0dot];
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
% lambda_min = -7.0e-13;
% lambda_max = -7.0e-14;
lambda_par = [lambda_min lambda_max];

param.input.lambda = lambda_par;

save(currentfile, 'param', '-append');

lambda_points = round((lambda_max-lambda_min)/delta_lambda0)

%% Chosing lambda_num

load(currentfile, 'param');

lambda_num = 8

param.input.lambda_num = lambda_num;

save(currentfile, 'param', '-append');

%% Counting photons

load(currentfile, 'param');

mp_ra = 0;
mp_dec = 0;

tic

('SELEZIONA I FOTONI DELLA SORGENTE')
%Using source_photon function to compute t of photons coming from a fixed direction
[t_fotoni_mjd, ra_deg, dec_deg] = source_photon(param.input.ra, param.input.dec, mp_ra, mp_dec, param.input.PEPOCH, param.input.ROI, param.input.Tstart, param.input.Tstop, param.input.E_min);
param.t_fotoni_mjd = t_fotoni_mjd;
param.ra_deg = ra_deg;
param.dec_deg = dec_deg;

save(currentfile, 'param', '-append');

toc
