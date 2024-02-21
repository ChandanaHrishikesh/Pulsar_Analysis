%% Creating runs

addpath(genpath('/data/virgo/Matlab/Snag/'))

currentfile = 'dB_current_RUN.mat';

load(currentfile, 'foldername', 'param');

d = datevec(datetime);
RUNfoldername = sprintf("%02d%02d%02dT%02d%02d", d(1)-2000, d(2), d(3), d(4), d(5));
path = fullfile(foldername, RUNfoldername)
mkdir(path)

Tp = param.input.Tp;
ROI = param.input.ROI;
n_pos = param.input.n_pos;
dt = param.input.dt;

epsilon = 0.54*2*dt*(7.0/Tp);
precision_sky = -dB_oomp(epsilon) + 3;
epsilon = round(epsilon, precision_sky);
ra_deg = round(param.ra_deg, precision_sky);
dec_deg = round(param.dec_deg, precision_sky);

if n_pos == 0
    n_pos = 2*ceil(ROI/epsilon) + 1
end
param.n_pos = n_pos;
save(currentfile, 'param', '-append');

str_skyS = strings(n_pos);

skySname = sprintf("SKYS.sh")
fID_skyS = fopen(fullfile(path, skySname), 'w');

tic

for k = 1:n_pos
    RA_grid = round(ra_deg + (k - (n_pos+1)/2)*epsilon, precision_sky);
    for kk = 1:n_pos
        DEC_grid = round(dec_deg + (kk - (n_pos+1)/2)*epsilon, precision_sky);

        filename_sky = sprintf("sky_%d_%d.mat", k, kk);
        str_skyS(k,kk) = filename_sky;

        R2 = (RA_grid - ra_deg)^2 + (DEC_grid - dec_deg)^2;

        if R2 <= ROI^2
            save(fullfile(path, filename_sky), 'RA_grid', 'DEC_grid', '-v7.3');
            corrname = sprintf("corr_%d_%d.m", k, kk)
            fID_corr = fopen(fullfile(path, corrname), 'w');
            fprintf(fID_corr, 'addpath(genpath(''/data/virgo/Work/All_new1/''))\n');
            fprintf(fID_corr, 'addpath(genpath(''/data/virgo/Matlab/Snag/''))\n');
            fprintf(fID_corr, 'dB_sky_corr(''%s'', ''%s'')\n', currentfile, fullfile(path, filename_sky));
            fclose(fID_corr);
            logname = sprintf("out_corr_%d_%d.log", k, kk);
            command = sprintf('nohup matlab -nodisplay -nosplash -r "run(''%s''); exit" > %s 2>&1 </dev/null &\n', corrname, logname);
            fprintf(fID_skyS, command);
        end
    end
end

fclose(fID_skyS);

save(currentfile, 'path', 'str_skyS', '-append');

toc
