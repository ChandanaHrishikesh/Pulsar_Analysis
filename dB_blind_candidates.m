function dB_blind_candidates(currentfile, skyfile, partfile)

load(currentfile, 'param')
load(skyfile, 't_corr', 'RA_grid', 'DEC_grid', 'delta_lambda', 'precision')
load(partfile, 'lambda_start')

Tp = param.input.Tp;
dt = param.input.dt;
lambda_num = param.input.lambda_num;
lambda_max = param.input.lambda(2);
CRThr = param.input.CRTH;

dts = t_corr - t_corr(1);

lambda = zeros(1,lambda_num);
CandS = zeros(6,ceil(Tp*86400/dt)*lambda_num);
Ncand = 0;

tic

for j = 1:lambda_num
    j
    lambda_new = round(lambda_start + (j-1)*delta_lambda, precision);
    if lambda_new <= lambda_max + delta_lambda/2
        lambda(j) = lambda_new;
        sdt = lambda(j)*dts.^2;
        t = t_corr + sdt;
        % ('SPETTRO-AUTOCOR')
        [~, sp_re, ~] = p_autcor_blind(t, Tp, dt); % sp_re e' meta' intervallo (sp_re is half range)

        y = y_gd(sp_re);
        f = x_gd(sp_re);
        df = dx_gd(sp_re);

        nn = round(0.01/df);
        y = y(nn:length(y));
        f = f(nn:length(f));

        Stat = robstat(y, 0.1);
        [Index] = candidate_selection_thr((y-Stat(1))/Stat(2), CRThr);

        s = sum(Index);
        CandS(1,Ncand+1:Ncand+s) = f(Index); % frequency
        CandS(2,Ncand+1:Ncand+s) = RA_grid;
        CandS(3,Ncand+1:Ncand+s) = DEC_grid;
        CandS(4,Ncand+1:Ncand+s) = 2*f(Index)*lambda(j); % Spin_dowmn
        CandS(5,Ncand+1:Ncand+s) = y(Index);
        CandS(6,Ncand+1:Ncand+s) = (y(Index)-Stat(1))/Stat(2);
        Ncand = Ncand + s;
    end
end

CandS = CandS(:,1:Ncand);

save(partfile, 'lambda', 'CandS', '-append');

clear lambda
clear CandS

toc
