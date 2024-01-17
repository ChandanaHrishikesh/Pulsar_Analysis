function [t_fotoni_mjd ra dec]= source_photon(ra,dec,mp_ra,mp_dec,PEPOCH,ROI,T_start,T_stop,E_min)


format long;
%KMAU= 1.49597870691000015e+08;  % conversione Km-AU
c_light= 299792.458;            % velocita' della luce in Km/s

rif_gps=mjd2gps(51910);               % mjd  di rifermento convertito in gps


%************* POSIZIONE DELLA SORGENTE********************************

ra=hour2deg(ra);


if dec(1)<0
   dec=dec(1)-dec(2)/60-dec(3)/3600;
else
    dec=dec(1)+dec(2)/60+dec(3)/3600;
end

%******** MOTO PROPRIO DELLA SORGENTE**************

t_medio=(T_stop - T_start)/2;  % t in mjd
t_mp=T_start+t_medio;          % tempo in mjd per il moto proprio
deltaT=mjd2gps(t_mp-PEPOCH);   % dt in gps

 
ra=  ra + mp_ra*(deltaT)/(86400*365.25*1000*3600);
dec= dec+ mp_dec*(deltaT)/(86400*365.25*1000*3600);
% 
% ************************************************************************

ra0=[ra-ROI ra+ROI] ;       % intorno di RA 
dec0=[dec-ROI dec+ROI];     % intorno di DEC

n=astro2rect([ra dec],0);   % vettore che punta la sorgente

%********************************************************

win=2;


%************************* MODIFICA 5/12/2012 *****************************
%**************************************************************************
% LE LISTE

% P6V11
% list_photon_p6v11_mat.txt CONTIENE TUTTI PARAMENTRI DEI DATI p6v11
% list_photon_reduced_p6v11_mat.txt CONTIENE SOLO I PARAMENTRI UTILIZZATI NELL'ANALISI DEI DATI p6v11

% P7V6
% list_photon_reduced_p7v6_mat.txt   CONTIENE SOLO I PARAMENTRI UTILIZZATI NELL'ANALISI DEI DATI p7v6

%**************************************************************************

%fid=fopen('/storage/users/tringali/ANALISI_PULSAR/LISTE/list_photon_p6v11_mat.txt','r'); % DATI DA P6V11 completi


%********************* LISTE RIDOTTE  ***********************************

% DATI DA P6V11
fid=fopen('/storage/users/tringali/FERMI_FILES/LISTE/list_photon_reduced_p6v11_mat.txt','r'); 

% DATI DA P7V6
%fid=fopen('/storage/users/tringali/FERMI_FILES/LISTE/list_photon_reduced_p7v6_mat.txt','r');

%**************************************************************************
%**************************************************************************



if fid < 0
    disp([list_photon ' file could not be opened']);
    return;
end
nfiles=0;

while (feof(fid) ~= 1)
    nfiles=nfiles+1;
    file{nfiles}=fscanf(fid,'%s',1);
    %str=sprintf('  %s ',file{nfiles});
    %disp(str);    
  
end

fclose(fid);

%******************* TEMPI INIZILI E FINALI DI P6V11 E P7V6***************

%  questa tabella e' data dal programma /ANALISI/PULSAR/read_week_photon.mat
% il programma read_week_photon.mat tira fuoriun file.txt con il tempo iniziale e finale di ogni settimana

% TEMPI DI P6V11 
 [t_ini t_fin]=textread('/storage/users/tringali/ANALISI_PULSAR/tempi_fotoni.txt','%f\t %f\t'); 

% TEMPI DI P7V6 
%[t_ini t_fin ]=textread('/storage/users/tringali/ANALISI_PULSAR/tempi_fotoni_p7v6.txt','%f\t %f\t'); 

%**************************************************************************


[aa,ini]=min(abs(t_ini-T_start));

[cc,fine]=min(abs(t_fin-T_stop));



for i = ini : fine
    fil=file{i};      
    loafil=['load ''' fil '''']; 
    eval(loafil)

    %t_fotoni=pdata{10};
    t_fotoni=pdata_red{5};
    
    %E_photon=pdata{1};      % energia fotoni MeV
    E_photon=pdata_red{1};
    
    a=pdata_red{2};       % ascensione della sorgente

    d=pdata_red{3};       % declinazione della sorgente
    
    zenit_angle= pdata_red{4};   % angolo zenitale
    
    %****************SELEZIONE in MJD SUI TEMPI di ARRIVO DEI FOTONI della PULSAR **************** 
    
    t_fotoni_mjd=gps2mjd(rif_gps+t_fotoni);    % t in MJD
    
%     jj=find(t_fotoni_mjd > T_start  & t_fotoni_mjd < T_stop);
%     
%     t_fotoni_mjd= t_fotoni_mjd(jj);
%     
%     E_photon=E_photon(jj);      % energia fotoni MeV
%     
%     a=a(jj);           % ascensione della sorgente
%     
%     d=d(jj);
    
    %************SELEZIONE IN DIREZIONE**********************************
    
    [iii]=fm_selezionesky_blind(E_photon,a,d,ra0,dec0,win);
       
    t_fotoni_mjd=t_fotoni_mjd(iii);  % t selezionati in (ra,dec) in MJD
   
    E_photon=E_photon(iii);
    
    zenit_angle=zenit_angle(iii);
 
    %************SELEZIONE SULL'ANGOLO ZENITALE****************************
    
    zz= find (zenit_angle < 105);
    
    t_fotoni_mjd=t_fotoni_mjd(zz);  % t selezionati in (ra,dec) in MJD
   
    E_photon=E_photon(zz);
    
    %******************SELEZIONE SULL'ENERGIA*******************************
    
    ee= find(E_photon > E_min);
    
    t_fotoni_mjd=t_fotoni_mjd(ee);
    
    %********************CORREZIONE LEAPSECOND ****************************
    
    
    t_fotoni_mjd=mjd2tai(t_fotoni_mjd);  % t+leapsecond in MJD
    
    t_fotoni_nocor= mjd2gps(t_fotoni_mjd); 
    
    
    %********* VETTORE DEI TEMPI CORRETTI DELLA SORGENTE*******************
    
    t_fotoni_mjd=t_fotoni_mjd';
    
      if i > 1
            t_fotoni_mjd=[t_ph t_fotoni_mjd];
      end
      
    t_ph=t_fotoni_mjd;
        
        
end



%********************CORREZIONE DELAY EINSTEIN ****************************
%('D-EINSTEIN')

t_deinstein=tdt2tdb(t_fotoni_mjd);

t_f_gps=mjd2gps(t_fotoni_mjd)+t_deinstein;

t_fotoni_mjd= gps2mjd(t_f_gps);