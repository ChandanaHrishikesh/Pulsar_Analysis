function dB_blind_candidates(file)

load(file, 't_corr', 'RA_grid', 'DEC_grid', 'Tp', 'dt', 'lambda_start', 'lambda_num', 'delta_lambda', 'lambda_max', 'CRThr', 'precision')

dts = t_corr - t_corr(1);

lambda = zeros(1,lambda_num);
CandS = [];

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

        %Ncand = length(Index);
        Cand(1,:) = f(Index); % frequency
        Cand(2,:) = RA_grid;
        Cand(3,:) = DEC_grid;
        Cand(4,:) = 2*f(Index)*lambda(j); % Spin_dowmn
        Cand(5,:) = y(Index);
        Cand(6,:) = (y(Index)-Stat(1))/Stat(2);
        CandS = [CandS Cand];
        clear Cand
    end
end

save(file, 'lambda', 'CandS', '-v7.3');

clear lambda
clear CandS

toc
