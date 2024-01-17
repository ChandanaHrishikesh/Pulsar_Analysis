function [p_dec_terra p_dec_sc t]=pos_detector_blind(t)



%function [p_dec_terra p_dec_sc t]=pos_detector(t)
%POS_DETECTOR  computes the position of the detector
%
%     doptab    table containing the Doppler data (depends on antenna)
%                (a matrix (8,n) 
%
%     t         time array(in mjd days)
%
% Converts UT in TT to obtain correct results

% Version 2.0 - December 2009
% Part of Snag toolbox - Signal and Noise for Gravitational Antennas
% by Sergio Frasca - sergio.frasca@roma1.infn.it
% Department of Physics - Universit� "Sapienza" - Rome



% ******************************MODIFICA ******************************
%   
%   SE VAR=1 CALCOLA LE SPLINE DELLA TERRA-SSB                        *
%   IN USCITA HO (X,Y,Z) IN AU                                        *
%   [p_dec t]=pos_detector([54600:10/86400:54800],1);                 *
% 
%   SE VAR !=1 CALCOLA LE SPLINE SPACECRAFT-TERRA
%   IN USCITA HO (X,Y,Z) IN METRI
%   [p_dec t]=pos_detector([54682.75054:10/86400:54762],2);
%
%
%    ALL'INIZIO LEGGE I TEMPI DI ARRIVO DEI FOTONI, TEMPI SUI QUALI VERRANNO
%    FATTE LE SPLINE.
% *********************************************************************


%deg_to_rad=pi/180;
rif_gps=mjd2gps(51910); 
KMAU= 1.49597870691000015e+08;  % conversione Km-AU

% *********************************************************************


% fid=fopen('/storage/users/tringali/Doppler/list_photon_mat_modif.txt','r');
% 
%     if fid < 0
%     disp([list_photon ' file could not be opened']);
%     return;
% end
% nfiles=0;
% 
% while (feof(fid) ~= 1)
%     nfiles=nfiles+1;
%     file{nfiles}=fscanf(fid,'%s',1);
%     str=sprintf('  %s ',file{nfiles});
%     disp(str);    
%   
% end
% 
% 
% fclose (fid)
% 
% 
% for i = 1 : nfiles-1
%     fil=file{i};      
%     loafil=['load ''' fil '''']; %
%     eval(loafil)
%     
%     t_fotoni=pdata{10};                        %tempi dei fotoni in GPS
%     t_fotoni_mjd = gps2mjd(rif_gps+t_fotoni);     % tempi dei fotoni in MJD
%  
%     t= t_fotoni_mjd;
%     if i > 1
%         t=[tpast;t];
%     end
%     tpast=t;
% end
% 
%     t=t';
% whos t;

%t=mjd2tt(t);  % arrya cge passo

% *********************************************************************

%p_dec=zeros(3,length(t));



%**************************** TERRA-SSB *****************************************

% terra_SSB=fopen('terra_ssb.txt','w');

[t_terra t_jps x y z ]=textread('/storage/users/tringali/pss_astro/table.dat','%f\t %f\t %f\t %f\t %f\t');
    
x= x*KMAU;
y= y*KMAU;
z= z*KMAU;
    
doptab_terra=[t_terra' ; x'; y'; z']; % ho messo nella matrice doptab i tempi della terra e le posizioni

%[n1,n2]=size(doptab_terra);


% tmin=min(t1)-0.5;
% tmax=max(t1)+0.5;
% doptab=reduce_doptab(doptab,tmin,tmax);



% p_dec(1:l,1)=spline(doptab(1:l,1),doptab(1:l,2),t1); % posizione x
% p_dec(1:l,2)=spline(doptab(1:l,1),doptab(1:l,3),t1); % posizione y
% p_dec(1:l,3)=spline(doptab(1:l,1),doptab(1:l,4),t1); % posizione z


    
p_dec_terra(1,:)=spline(doptab_terra(1,:),doptab_terra(2,:),t);    % posizione x-spline in AU
p_dec_terra(2,:)=spline(doptab_terra(1,:),doptab_terra(3,:),t);    % posiZione y-spline in AU
p_dec_terra(3,:)=spline(doptab_terra(1,:),doptab_terra(4,:),t);    % posizione z-spline in AU


% p_dec(1,:)=p_dec(1,:)*KMAU;    % posizione x-spline in km
% p_dec(2,:)=p_dec(2,:)*KMAU;    % posizione y-spline in km
% p_dec(3,:)=p_dec(3,:)*KMAU;    % posizione z-spline in km

    


% fprintf(terra_SSB,'%15.10f\t %10.10f\t %10.10f\t %10.10f\t \n',[t;p_dec]);

%     figure;hold on;
%     plot(doptab_terra(1,:),doptab_terra(2,:),'ro');
%     plot(t,p_dec_terra(1,:),'*');
%     title('TERRA')

% **************************** SPACECRAFT-TERRA **************************

%************************* MODIFICA 5/12/2012 ************************
% LE LISTE:

% list_spacecraft_mat.txt CONTIENE TUTTI PARAMENTRI
% list_spacecraft_reduced_mat.txt CONTIENE SOLO I PARAMENTRI UTILIZZATI NELL'ANALISI


('LEGGO IL FILE DI FERMI')
   
%fid=fopen('/storage/users/tringali/ANALISI_PULSAR/LISTE/list_spacecraft_mat.txt','r');

fid=fopen('/storage/users/tringali/FERMI_FILES/LISTE/list_spacecraft_p6v11_reduced_mat.txt','r'); %  p6v11


%***************  BLOCCO DI DATI P7V6*****************************

%fid=fopen('/storage/users/tringali/FERMI_FILES/LISTE/list_spacecraft_p7v6_reduced_mat.txt','r');   


%*************************************************************************

%*************************************************************************

% USO QUESTO fid QUANDO FACCIO GIRARE IL PROGRAMMA spettro_spindown.m
%fid=fopen('/storage/users/tringali/spindown/fermi.txt','r');


%     if fid < 0
%         disp([list_sc ' file could not be opened']);
%         return;
%     end

nfiles=0;

while (feof(fid) ~= 1)
    
    nfiles=nfiles+1;
    file{nfiles}=fscanf(fid,'%s',1); % leggi da file il nome del file.mat
%     str=sprintf('  %s ',file{nfiles});% stampa il nome del file-settimana.mat
%     disp(str);    
%   
end


fclose(fid);

%sc_terra=fopen('sc_terra.txt','w');

for i = 1 : nfiles-1
    
    fil=file{i};       %prende il nome del filei-simo della lista caricata prima
    loafil=['load ''' fil '''']; %
    eval(loafil)
    
    %t_sc=pdata{1};                        %tempi dello spacecraft in GPS
    t_sc=pdata_red{1};
    
    t_sc_mjd = gps2mjd(rif_gps+t_sc);     % tempi dello spacecraft in MJD
    
    %pos=pdata{3};         %posizioni dello spacecraft in metri
    pos=pdata_red{2}; 
    
    
    
    pos_x = pos(:,1)*1e-03;
    pos_y = pos(:,2)*1e-03;
    pos_z = pos(:,3)*1e-03;
    
    rock_angle= pdata_red{3};   % rock-angle dello spacecraft
    lat_config = pdata_red{4};
    data_qual = pdata_red{5};
    
    %   ********************** SELEZIONE SUL ROCK-ANCGLE********************
    
    j=find(rock_angle<52);
    
    t_sc_mjd=t_sc_mjd(j);
    pos_x=pos_x(j);
    pos_y=pos_y(j);
    pos_z=pos_z(j);
    lat_config=lat_config(j);
    data_qual=data_qual(j);
    
    %**********************************************************************
    
    %*********************** LAT_CONFIG E DATA_QUAL************************
    
    jj=find(lat_config ==1 & data_qual==1);

    t_sc_mjd=t_sc_mjd(jj);
    pos_x=pos_x(jj);
    pos_y=pos_y(jj);
    pos_z=pos_z(jj);
 
    %**********************************************************************
        
	doptab_sc=[t_sc_mjd' ; pos_x'; pos_y'; pos_z'];
    
	if i > 1
        doptab_sc=[doppast doptab_sc];
    end
    
    doppast=doptab_sc;
end
    
    
p_dec_sc(1,:)=spline(doptab_sc(1,:),doptab_sc(2,:),t);    % posizione x-spline in m
p_dec_sc(2,:)=spline(doptab_sc(1,:),doptab_sc(3,:),t);    % posiZione y-spline in m
p_dec_sc(3,:)=spline(doptab_sc(1,:),doptab_sc(4,:),t);    % posizione z-spline in m
    

%     p_dec(1,:)=p_dec(1,:)*1e-03;    % posizione x-spline in km
%     p_dec(2,:)=p_dec(2,:)*1e-03;    % posizione y-spline in km
%     p_dec(3,:)=p_dec(3,:)*1e-03;    % posizione z-spline in km
    
   % fprintf(sc_terra,'%15.10f\t %15.10f\t %15.10f\t %15.10f\t\n',[t;p_dec]);

   

% 
%     figure;hold on;
%     plot(doptab_sc(1,:),doptab_sc(2,:),'ro');
%     plot(t,p_dec_sc(1,:),'*');
%     title('SPACECRAFT')

    


