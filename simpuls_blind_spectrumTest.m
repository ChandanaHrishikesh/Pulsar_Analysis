function [ampiez_cand_max freq_cand_max f0dot_cand_max lambda_fin_max RA_cand_max DEC_cand_max]=simpuls_blind_spectrumTest(ra,dec,mp_ra,mp_dec,ROI,T_start,T_stop,E_min,PEPOCH,spin_par,n_pos,Tp,dt,lambda_par,blind)
%function [cl_grid x_grid ampiez_cand freq_cand f0dot_cand lambda_fin RA_cand DEC_cand]=simpuls_blind(ra,dec,mp_ra,mp_dec,ROI,T_start,T_stop,E_min,PEPOCH,spin_par,n_pos,n_f0dot,Tp,dt)
%function [ampiez_cand_max freq_cand_max f0dot_cand_max lambda_fin_max RA_cand_max DEC_cand_max]=simpuls_blind_spectrumTest(ra,dec,mp_ra,mp_dec,ROI,T_start,T_stop,E_min,PEPOCH,spin_par,n_pos,Delta_f1,Tp,dt,numrun)

% SORGENTE J0633+0632:
% f0dot=-8.999148066954602476  10?^-13
% freq=3.3625291587551729515
%[ampiez_cand freq_cand f0dot_cand lambda_fin RA_cand DEC_cand]=simpuls_blind_spectrum('06:33:44.2159045',[06 32 34.94512],0,0,0.6,54682.655439529444266,59624.61388888889,550,54945,[0 -9.0e-13 0],6,10,7,1/64)

% PARAMETRI in INPUT (Input parameters)
% 
% ra ascensione retta della sorgente in ore(ra right ascension of the source in hours)
% dec declinazione della sorgente in ore ma deve esser data come un
% vettore(dec declination of source in hours but must be given as a vector)
%     
% 
% ROI raggio,intorno in gradi usato per (ra,dec) (ROI radius, around in degrees used for (ra, dec))
% 
% T_start tempo di inizio dei dati relativi alla pulsar (T_start start time of pulsar data)
% T_stop  tempo di fine dei dati relativi alla pulsar(T_stop end time of pulsar data)
% 
% E_min energia minima (E_min minimum energy)
% 
% PEPOCH epoca di riferimento della pulsar in MJD (PEPOCH reference epoch of the pulsar in MJD)
% 
% spin_par  e' un vettore che contiene i parametri di spin-down [f0 f1 f2]
% (spin_par is a vector containing the spin-down parameters [f0 f1 f2])


tic

% if exist('RunNumber.dat','file')==2     % Load the previous run number if it exists, or initialise it to 0
%     load('RunNumber.dat','RunNumber');
% else
%     RunNumber=0;
% end
% 
% RunNumber=RunNumber+1;   %Increment the run number
% save('RunNumber.dat','RunNumber');
% 
% RunDateTime=datetime('now');
format long;

%N=round(Tp*86400/dt);


SP_CAND=[];
POS=[];
RA_MAX=[];
DEC_MAX=[];          

RA_cand_max=0; %first time when you have the value 
DEC_cand_max=0;




('SELEZIONA I FOTONI DELLA SORGENTE')
%Using source_photon function
[t_fotoni_mjd ra_deg dec_deg]= source_photon(ra,dec,mp_ra,mp_dec,PEPOCH,ROI,T_start,T_stop,E_min); % t dei fotoni che vengono da una fissata direzione (t of photons coming from a fixed direction)

h=0;

save PIPPO
pause

%epsilon= 0.54/32;           % passo della griglia (grid pitch)

epsilon= 0.54*2*dt*(7.0/Tp);
%for i = 1:numrun;
ra_in=ra_deg;
dec_in=dec_deg;
%ra_deg=ra_deg+(epsilon*normrnd(0,1));  %central point
%dec_deg=dec_deg+(epsilon*normrnd(0,1));

alpha = ra_deg-(n_pos/2)*epsilon;   %ascensione retta in gradi da cui parto per costruire la griglia (right ascension in degrees from which I start to build the grid)

delta = dec_deg-(n_pos/2)*epsilon;  % declinazione in gradi da cui parto per costruire la griglia (declination in degrees from which I start to build the grid)


('COSTRUZIONE DELLA GRIGLIA IN POSIZIONE E SPINDOWN')

%* CICLO su RA, DEC e PARAMETRI di SPIN-DOWN'(CYCLE on RA, DEC and SPIN-DOWN' PARAMETERS)


k_min = 0;
k_max = n_pos;

kk_min = 0;
kk_max = n_pos;
fid=fopen('Check2.txt','a+');
%filename=sprintf('selected_30d_30_0633_0632.mat');
%filename=sprintf('selected_30d_30_0.003938_0633_0632.mat');
%filename=sprintf('selected_new_%dd_%d_%s%s%s%s_%02d%02d.mat',Tp,n_pos,ra(1),ra(2),ra(4),ra(5),dec(1),dec(2));
%file_image=sprintf('spectrum_new_%dd_%d_%s%s%s%s_%02d%02d.jpg',Tp,n_pos,ra(1),ra(2),ra(4),ra(5),dec(1),dec(2));
%Tobs
filename=sprintf('Changed__%dd_%d_%f_%s%s%s%s_%02d%02d.mat',Tp,n_pos,T_stop,ra(1),ra(2),ra(4),ra(5),dec(1),dec(2));
file_image=sprintf('spectrum_chan_%dd_%d_%f_%s%s%s%s_%02d%02d.jpg',Tp,n_pos,T_stop,ra(1),ra(2),ra(4),ra(5),dec(1),dec(2));
hh=0;
n_grid=astro2rect([ra_deg dec_deg],0);
[t_corr]= doppler_correction_blind (t_fotoni_mjd, n_grid);
t=t_corr;                      % tempi corretti per il Doppler in GPS e per lp SHAPIRO (times corrected for Doppler in GPS and for SHAPIRO)
T_obs = t(length(t)) - t(1);
%lambda_input= f1dot/(2*f0);           % Adding line for lambda_input
delta_f0dot = 1/(2* T_obs*Tp*86400);
delta_lambda = (delta_f0dot*2*dt)/2; % c'e' 2*dt perche' la frequenza massima e' 1/64 ma quella fisica, che ci interessa e' 1/32 (there is 2*dt because the maximum frequency is 1/64 but the physical one, which interests us, is 1/32)
%n_maria=abs(lambda_input)/delta_lambda;
fprintf(fid,'Filename=%s\n',filename);
fprintf(fid,'ra_in=%d\n dec_in=%d\n ra_deg=%d\n dec_deg=%d\n',ra_in, dec_in, ra_deg, dec_deg);% Central point and Inputs
%fprintf(fid,'n_maria=%d\n delta_lambda=%d\n',n_maria, delta_lambda);

%save PIPPO2
for k = k_min:k_max % sky loop
    
    RA_grid= alpha + k*epsilon;
   

    for kk = kk_min:kk_max 
        
        h=h+1;
        
        DEC_grid= delta + kk*epsilon;
        
        
        n_grid=astro2rect([RA_grid DEC_grid],0);
        
        
        [t_corr]= doppler_correction_blind (t_fotoni_mjd, n_grid);
 
        t=t_corr;                      % tempi corretti per il Doppler in GPS e per lp SHAPIRO (times corrected for Doppler in GPS and for SHAPIRO)

        T_obs = t(length(t)) - t(1);   % tempo di osservazione in secondi (observation time in seconds)
        
        % fprintf(fid,'T(%d,%d)=%f, ',k,kk,T_obs);


('CORREZIONE SPIN-DOWN')
% ****************** CORREZIONE DI SPIN-DOWN *************************
 
        tempo_primo_fotone=t(1); 

        t_pulsar=mjd2gps(PEPOCH);   % tempo di riferimento della pulsar in GPS (reference time of the pulsar in GPS)

        dts=t-tempo_primo_fotone;

        t_orig= t-tempo_primo_fotone;

%******* riscrivo i parametri rispetto al primo tempo di arrivo del
%fotone(I rewrite the parameters with respect to the first arrival time of the photon)
    
        if blind==0 % for known source parameters
        %if (spin_par(3) ~= 0)   
            f0=spin_par(1);
            f1dot=spin_par(2);
            f2dot=spin_par(3);
            f0= f0+f1dot*(tempo_primo_fotone - t_pulsar) + f2dot *(tempo_primo_fotone - t_pulsar)^2/2;
            f1dot=f1dot+f2dot*(tempo_primo_fotone - t_pulsar);
            sdt=f1dot*dts.^2/(2*f0)+f2dot*dts.^3/(6*f0);   % correzione al primo ordine (first order correction) 
            t=t+sdt;           
            delta_f0dot = 1/(2* T_obs*Tp*86400);
            delta_lambda = (delta_f0dot*2*dt)/2;
                       
        else % around a frequency and spin-down range defined by lambda_par
            lambda_min=lambda_par(1);
            lambda_max=lambda_par(2);
            delta_f0dot = 1/(2* T_obs*Tp*86400);
            delta_lambda = (delta_f0dot*2*dt)/2; % c'e' 2*dt perche' la frequenza massima e' 1/64 ma quella fisica, che ci interessa e' 1/32 (there is 2*dt because the maximum frequency is 1/64 but the physical one, which interests us, is 1/32)
            lambda_num=round((lambda_max-lambda_min)/delta_lambda);
           
            for j=1:lambda_num
                hh=hh+1;
                %lambda(j) = -470*delta_lambda-(j-1)*delta_lambda;
                %lambda(j) = -((n_maria+j-2))*delta_lambda;     % written in general form
                lambda(j) = lambda_min + (j-1)*delta_lambda;
                sdt = lambda(j)*dts.^2;
                t=t_corr+sdt;
                ('SPETTRO-AUTOCOR')
                [autc sp_re sp_im]=p_autcor_blind(t,Tp,dt); % sp_re e' meta' intervallo (sp_re is half range)

                y=y_gd(sp_re);
                f=x_gd(sp_re);
                df=dx_gd(sp_re);


                nn=round(0.01/df);
                y=y(nn:length(y));
                f=f(nn:length(f));

                [sp_cand(j) pos(j)]= max(y); % to be changed by a function to select all candidates above a thresh.

                % y_max=sp_cand(j);
                f_spcand=f(pos(j));
                %f0dot_sp= 2*f_spcand*lambda(j);

                np= 24;
                per=1/f_spcand;
                dx=per/np;
                x=dx/2:dx:per;
                x(length(x)+1)=x(length(x))+dx;

                [lc]=lightcurve_blind(t,per,np,x);
                cl_grid(hh,:)= lc;
                x_grid(hh,:)= x;
                %x_grid(hh,length(x)+1)=x(length(x))+dx;
            end

             [val1 pos1]= max(sp_cand);  %max dei valori massimi trovati per i diversi valori di lambda(max of the maximum values ​​found for the different lambda values)
 
             f0_cand=( pos(pos1)+nn-2 )*df;      % freq corrispondente (corresponding frequency)
             lambda_cand= lambda(pos1);  %lambda corrispondete (corresponding lambda)
             f0_dot_cand = 2*f0_cand*lambda_cand;   % f02dot del candidato
        end
        
                        
% %******COSTRUISCONO DEI VETTORI SU CUI MEMORIZZO I VALORI DEI PARAMETRI (This builds vectors on which I store the values of the parameters)****
         ampiez_cand(h) = val1;
         freq_cand(h) = f0_cand;
         f0dot_cand(h) = f0_dot_cand;
         lambda_fin(h) = lambda_cand;
         RA_cand(h) = RA_grid;
         DEC_cand(h) = DEC_grid;
                 
%*************************************************************************  
       % kk=kk+1;  This line is not required
    end
    
    fprintf(fid,'\n');
end
fprintf(fid,'df=%d\n',df);
[val_max pos_max]= max(ampiez_cand);

ampiez_cand_max = ampiez_cand(pos_max);
freq_cand_max=  freq_cand(pos_max);
f0dot_cand_max= f0dot_cand(pos_max);
lambda_fin_max = lambda_fin(pos_max);
RA_cand_max = RA_cand(pos_max);
DEC_cand_max = DEC_cand(pos_max);

k_cand_max = round((RA_cand_max - alpha)/epsilon);
kk_cand_max = round((DEC_cand_max - delta)/epsilon);
j_fin_max=round(1-lambda_fin_max/delta_lambda-(n_maria-1));

%lambda_check=(lambda_input-lambda_fin_max)/delta_lambda;
%fprintf(fid,'lambda_check=%d\n',lambda_check);
%fprintf(fid,'else delta_lambda=%d\n',delta_lambda);
save(filename,'ampiez_cand_max','freq_cand_max','f0dot_cand_max','lambda_fin_max','RA_cand_max','DEC_cand_max','k_cand_max','kk_cand_max','j_fin_max','lambda_par','delta_lambda');
fprintf(fid,'amp_cand_max=%d\n, freq_cand_max=%d\n, f0dot_cand_max=%d\n, lambda_fin_max=%d\n, RA_cand_max=%d\n, DEC_cand_max=%d\n, k_cand_max=%d\n, kk_cand_max=%d\n', ampiez_cand_max,freq_cand_max, f0dot_cand_max, lambda_fin_max, RA_cand_max, DEC_cand_max, k_cand_max, kk_cand_max);
fprintf(fid,'k=%d\n kk=%d\n', k, kk);% Best candidate point



('PLOTTING THE SPECTRUM')


n_grid=astro2rect([RA_cand_max DEC_cand_max],0);


[t_corr]= doppler_correction_blind (t_fotoni_mjd, n_grid);

t=t_corr;                      % tempi corretti per il Doppler in GPS e per lp SHAPIRO (times corrected for Doppler in GPS and for SHAPIRO)

T_obs = t(length(t)) - t(1);   % tempo di osservazione in secondi (observation time in seconds)


('CORREZIONE SPIN-DOWN')
% ****************** CORREZIONE DI SPIN-DOWN *************************

tempo_primo_fotone=t(1);

t_pulsar=mjd2gps(PEPOCH);   % tempo di riferimento della pulsar in GPS (reference time of the pulsar in GPS)

dts=t-tempo_primo_fotone;

t_orig= t-tempo_primo_fotone;

%******* riscrivo i parametri rispetto al primo tempo di arrivo del
%fotone(I rewrite the parameters with respect to the first arrival time of the photon)


sdt = lambda_fin_max*dts.^2;
t=t_corr+sdt;

('SPETTRO-AUTOCOR')
[autc sp_re sp_im]=p_autcor_blind(t,Tp,dt); % sp_re e' meta' intervallo (sp_re is half range)

y=y_gd(sp_re);
f=x_gd(sp_re);
df=dx_gd(sp_re);
%  fprintf(fid,'df=%d\n',df);

nn=round(0.01/df);
y=y(nn:length(y));
f=f(nn:length(f));
Median=median(y)
[Std,Mean]=std(y)
CR=(ampiez_cand_max - Median)/sqrt(Std^2+(Mean-Median)^2)
fprintf(fid,'CR=%s\n',CR);
figure,plot(f,y)
hold on
grid on
xlabel('Frequency')
ylabel('Amplitude')
title(sprintf('Spectrum for lambda = %g, j = %d, RA_grid = %d,\n k = %d, DEC_grid = %d, kk = %d' , lambda_fin_max, j_fin_max, RA_cand_max, k_cand_max, DEC_cand_max, kk_cand_max))
fclose(fid);
ax = gca;
exportgraphics(ax,file_image)
%end

% Data_table=table(RunNumber, RunDateTime, ra_in, dec_in, ra_deg, dec_deg, RA_cand_max, DEC_cand_max, freq_cand_max, ampiez_cand_max, 'VariableNames', {'RunNumber', 'Date_Time', 'Input_RA', 'Input_dec', 'Central_RA', 'Central_dec', 'Output_RA', 'Output_dec', 'Freq_max', 'Amp_max'});  % Table with information to be stored
% writetable(Data_table, 'Data_table.dat', 'Delimiter', '\t');
toc





  