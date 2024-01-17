function read_week_photon 


% il programma read_week_photon.mat tira fuori il tempo iniziale e finale di ogni settimana
% per ogni fotone



%fid=fopen('/storage/users/tringali/ANALISI_PULSAR/LISTE/list_photon_reduced_mat.txt','r'); % dati da p6v11 

fid=fopen('/storage/users/tringali/FERMI_FILES/LISTE/list_photon_reduced_p7v6_mat.txt','r'); % dati da p7v6

if fid < 0
    disp([listin ' file could not be opened'])
    return;
end
nfiles=0;

while (feof(fid) ~= 1)
    nfiles=nfiles+1;
    file{nfiles}=fscanf(fid,'%s',1);
    str=sprintf('  %s ',file{nfiles});
    disp(str);
end

fclose(fid);


%tempi_fotoni=fopen('tempi_fotoni.txt','w');

tempi_fotoni_p7v6=fopen('tempi_fotoni_p7v6.txt','w');

rif_gps=mjd2gps(51910);
t_ini=zeros(1,nfiles-1);
t_fin=t_ini;


for i = 1 : nfiles-1
   
    fil=file{i};      
    loafil=['load ''' fil '''']; 
    eval(loafil)
    
    
    %t_fotoni=pdata{10};   % tempi di arrivo dei fotoni in GPS
    t_fotoni=pdata_red{5}; % tempi di arrivo dei fotoni in GPS ma dal file.mat ridotto
    t_fotoni_mjd=gps2mjd(rif_gps+t_fotoni);
    
    
    t_ini(i)=t_fotoni_mjd(1);
    t_fin(i)=t_fotoni_mjd(length(t_fotoni_mjd));
    
    %fprintf(tempi_fotoni,'%15.10f\t %15.10f\n',t_ini(i),t_fin(i));
    
    fprintf(tempi_fotoni_p7v6,'%15.10f\t %15.10f\n',t_ini(i),t_fin(i));
end   