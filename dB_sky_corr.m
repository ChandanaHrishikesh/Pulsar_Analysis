function dB_sky_corr(currentfile, skyfile)

load(currentfile, 'param')
load(skyfile, 'RA_grid', 'DEC_grid')

Tp = param.input.Tp;
dt = param.input.dt;

t_fotoni_mjd = param.t_fotoni_mjd;

tic

('COMPUTING CORRECTIONS')

n_grid = astro2rect([RA_grid DEC_grid], 0);

[t_corr] = doppler_correction_blind(t_fotoni_mjd, n_grid); % tempi corretti per il Doppler in GPS e per lp SHAPIRO (times corrected for Doppler in GPS and for SHAPIRO)

T_obs = t_corr(length(t_corr)) - t_corr(1); % tempo di osservazione in secondi (observation time in seconds)
delta_f0dot = 1/(2*T_obs*Tp*86400);
delta_lambda = (delta_f0dot*2*dt)/2; % c'e' 2*dt perche' la frequenza massima e' 1/64 ma quella fisica, che ci interessa e' 1/32 (there is 2*dt because the maximum frequency is 1/64 but the physical one, which interests us, is 1/32)
precision = -dB_oomp(delta_lambda) + 3;
delta_lambda = round(delta_lambda, precision);

save(skyfile, 't_corr', 'T_obs', 'delta_lambda', 'precision', '-append');

toc
