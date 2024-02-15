%% Creating runs

addpath(genpath('/media/cta/EXTERNAL_USB/Chandana_Work/Snag/Snag'))

load('dB_current_RUN.mat', 'foldername', 'param');

d = datevec(datetime);
RUNfoldername = sprintf("%02d%02d%02dT%02d%02d", d(1)-2000, d(2), d(3), d(4), d(5));
path = fullfile(foldername, RUNfoldername)
mkdir(path)

Tp = param.input.Tp;
ROI = param.input.ROI;
n_pos = param.input.n_pos;
dt = param.input.dt;
lambda_min = param.input.lambda(1);
lambda_max = param.input.lambda(2);
CRThr = param.input.CRTH;
lambda_num = param.input.lambda_num;

t_fotoni_mjd = param.t_fotoni_mjd;

epsilon = 0.54*2*dt*(7.0/Tp);
precision_sky = -dB_oomp(epsilon) + 3;
epsilon = round(epsilon, precision_sky);
ra_deg = round(param.ra_deg, precision_sky);
dec_deg = round(param.dec_deg, precision_sky);

str_skyS = strings(n_pos);

h = 0;
r = 0;
lambda_part_max = 1;

tic

for k = 1:n_pos
    RA_grid = round(ra_deg + (k - (n_pos+1)/2)*epsilon, precision_sky);
    for kk = 1:n_pos
        DEC_grid = round(dec_deg + (kk - (n_pos+1)/2)*epsilon, precision_sky);
        
        R2 = (RA_grid - ra_deg)^2 + (DEC_grid - dec_deg)^2;

        if R2 <= ROI^2
            filename_sky = sprintf("sky_%d_%d.mat", k, kk)
            str_skyS(k,kk) = filename_sky;

            n_grid = astro2rect([RA_grid DEC_grid], 0);

            [t_corr] = doppler_correction_blind(t_fotoni_mjd, n_grid); % tempi corretti per il Doppler in GPS e per lp SHAPIRO (times corrected for Doppler in GPS and for SHAPIRO)

            T_obs = t_corr(length(t_corr)) - t_corr(1); % tempo di osservazione in secondi (observation time in seconds)
            delta_f0dot = 1/(2*T_obs*Tp*86400);
            delta_lambda = (delta_f0dot*2*dt)/2; % c'e' 2*dt perche' la frequenza massima e' 1/64 ma quella fisica, che ci interessa e' 1/32 (there is 2*dt because the maximum frequency is 1/64 but the physical one, which interests us, is 1/32)
            precision = -dB_oomp(delta_lambda) + 3;
            delta_lambda = round(delta_lambda, precision);

            r = r + 1;
            runSname = sprintf("runS_%d.sh", r)
            fID_runS = fopen(fullfile(path, runSname), 'w');

            lambda_part = 0;
            lambda_start = lambda_min;
            while lambda_start <= lambda_max + delta_lambda/2
                lambda_part = lambda_part + 1;
                h = h + 1;
                dataname_h = sprintf("data_%d_%d_%d_%d.mat", h, k, kk, lambda_part);
                str_dataS(h) = dataname_h;
                save(fullfile(path, dataname_h), 't_corr', 'RA_grid', 'DEC_grid', 'Tp', 'dt', 'lambda_start', 'lambda_num', 'delta_lambda', 'lambda_max', 'CRThr', 'precision', '-v7.3');
                runname_h = sprintf("run_%d.m", h)
                fID_run_h = fopen(fullfile(path, runname_h), 'w');
                fprintf(fID_run_h, 'addpath(genpath(''/media/cta/EXTERNAL_USB/Chandana_Work/Pulsar_analysis''))\n');
                fprintf(fID_run_h, 'addpath(genpath(''/media/cta/EXTERNAL_USB/Chandana_Work/Snag/Snag''))\n');
                fprintf(fID_run_h, 'dB_blind_candidates(''%s'')\n', dataname_h);
                fclose(fID_run_h);
                logname_h = sprintf("out_run_%d.log", h)
                command_h = sprintf('nohup matlab -nodisplay -nosplash -r "run(''%s''); exit" > %s 2>&1 &\n', runname_h, logname_h);
                fprintf(fID_runS, command_h);
                lambda_start = round(lambda_start + lambda_num*delta_lambda, precision);
            end
            save(fullfile(path, filename_sky), 'T_obs', 'delta_lambda', 'lambda_part', '-v7.3');
            lambda_part_max = max(lambda_part_max, lambda_part);
            fclose(fID_runS);
        end
    end
end

param.lambda_part_max = lambda_part_max;

save('dB_current_RUN.mat', 'foldername', 'RUNfoldername', 'path', 'str_skyS', 'str_dataS', 'param', '-v7.3');

toc
