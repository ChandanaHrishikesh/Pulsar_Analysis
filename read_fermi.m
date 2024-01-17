function file=read_fermi(listin)

% QUESTO PROGRAMMA LEGGE LA LISTA.TXT DEI FILE.FITS E CREA I FILE.MAT
% listin="hjghj.txt" % %filename%
fid=fopen(listin,'r');
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

for i = 1 : nfiles-1

    filename=file{i}
    [pathstr, name, ext] = fileparts(filename);
    pdata = fitsread(filename,'bintable');
   % save(name,'pdata_red'),i
    
    
    %********* MODIFICA 4/12/2012 *****************************************
    % 
    % DAI FILE.FITS SELEZIONIAMO SOLO I PARAMETRI CHE VENGONO UTILIZZATI
    %NELL'ANALISI
    
    %*********** FOTONI ***************************************************
    
    pdata_red{1}= pdata{1};   % energia dei fotoni
    pdata_red{2}= pdata{2};   % ascensione retta
    pdata_red{3}= pdata{3};   % declinazione
    pdata_red{4}= pdata{8};   % angolo zenitale
    pdata_red{5}= pdata{10};  % tempo  di arrivo dei fotoni
    %**********************************************************************
    
    %********* SPACECRAFT**************************************************
     
%     pdata_red{1}= pdata{1};   % tempo dello spacecraft
%     pdata_red{2}= pdata{3};   % posizione dello spacecraft
%     pdata_red{3}= pdata{19};  % rock-angle dello spacecraft
%     pdata_red{4}= pdata{21};  % lat_config dello spacecraft
%     pdata_red{5}= pdata{22};  %data_qual dello spacecraft
    %**********************************************************************
     
    save(name,'pdata_red'),i
end
