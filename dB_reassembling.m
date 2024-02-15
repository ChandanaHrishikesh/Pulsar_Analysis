%% Reassembling results

load('dB_current_RUN.mat', 'foldername', 'RUNfoldername', 'path', 'str_skyS', 'str_dataS', 'param');

n_pos = param.input.n_pos;
lambda_num = param.input.lambda_num;
lambda_part_max = param.lambda_part_max;

h = 0;
T_obsS = zeros(n_pos,n_pos);
delta_lambdaS = zeros(n_pos,n_pos);
lambda_partS = zeros(n_pos,n_pos);
lambdaS = zeros(n_pos*lambda_num,n_pos*lambda_part_max);
CAND = [];

tic

for k = 1:n_pos
    for kk = 1:n_pos
        filename_sky = sprintf(fullfile(path, str_skyS(k,kk)))
        if isfile(filename_sky)
            load(filename_sky);
            T_obsS(k,kk) = T_obs;
            delta_lambdaS(k,kk) = delta_lambda;
            lambda_partS(k,kk) = lambda_part;
            for n = 1:lambda_part
                h = h + 1;
                dataname_h = sprintf(fullfile(path, str_dataS(h)));
                load(dataname_h);
                lambdaS(1+(k-1)*length(lambda):k*length(lambda),(kk-1)*lambda_part_max+n) = lambda.';
                CAND = [CAND CandS];
            end
        end
    end
end

param.T_obsS = T_obsS;
param.delta_lambdaS = delta_lambdaS;
param.lambda_partS = lambda_partS;
param.lambdaS = lambdaS;

CANDfilename = sprintf('CAND.mat');

save(fullfile(path, CANDfilename), 'CAND', 'param', '-v7.3');

toc

REF = load(fullfile(foldername, 'CAND_REF.mat'));
load(fullfile(path, CANDfilename))

('DONE')
