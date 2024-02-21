%% Creating runs

addpath(genpath('/data/virgo/Work/All_new1/'))

currentfile = 'dB_current_RUN.mat';

load(currentfile, 'path', 'str_skyS', 'param');

Tp = param.input.Tp;
dt = param.input.dt;
lambda_min = param.input.lambda(1);
lambda_max = param.input.lambda(2);
CRThr = param.input.CRTH;
lambda_num = param.input.lambda_num;

n_pos = param.n_pos;

h = 0;
r = 0;
lambda_part_max = 1;

tic

for k = 1:n_pos
    for kk = 1:n_pos
        fullfilename_sky = sprintf(fullfile(path, str_skyS(k,kk)))
        if isfile(fullfilename_sky)
            load(fullfilename_sky, 'delta_lambda', 'precision');
            
            r = r + 1;
            runSname = sprintf("RUNS_%d.sh", r)
            fID_runS = fopen(fullfile(path, runSname), 'w');

            lambda_part = 0;
            lambda_start = lambda_min;
            while lambda_start <= lambda_max + delta_lambda/2
                lambda_part = lambda_part + 1;
                save(fullfilename_sky, 'lambda_part', '-append');
                index = [k, kk, lambda_part];
                h = h + 1;
                dataname = sprintf("data_%d.mat", h);
                str_dataS(h) = dataname;
                save(fullfile(path, dataname), 'index', 'lambda_start', '-v7.3');
                runname = sprintf("run_%d.m", h)
                fID_run = fopen(fullfile(path, runname), 'w');
                fprintf(fID_run, 'addpath(genpath(''/data/virgo/Work/All_new1''))\n');
                fprintf(fID_run, 'addpath(genpath(''/data/virgo/Work/All_new1/''))\n');
                fprintf(fID_run, 'dB_blind_candidates(''%s'', ''%s'', ''%s'')\n', currentfile, fullfilename_sky, dataname);
                fclose(fID_run);
                logname = sprintf("out_run_%d.log", h);
                command = sprintf('nohup matlab -nodisplay -nosplash -r "run(''%s''); exit" > %s 2>&1 </dev/null &\n', runname, logname);
                fprintf(fID_runS, command);
                lambda_start = round(lambda_start + lambda_num*delta_lambda, precision);
            end
            lambda_part_max = max(lambda_part_max, lambda_part);
            fclose(fID_runS);
        end
    end
end

param.lambda_part_max = lambda_part_max;

save(currentfile, 'str_dataS', 'param', '-append');

toc
