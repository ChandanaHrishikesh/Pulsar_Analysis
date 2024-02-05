function [iii A mua mud sda sdd a d]=fm_selezionesky(e,a,d,ra0,dec0,win)


%function [ii A mua mud sda sdd a d]=fm_selezionesky(fil,ra0,dec0,win)
%fm_selsky   selection of events in sky
%
% Creare la lista con dir /s/b > listphot.dat ed editare
%
%   ra0,dec0  [min max]
%   nsig      number of sigma

% IL PROGRAMMA E' STATO MODIFICATO PER POTER ESSER UTILIZZATO 
% NELLA CORREZIONE DOPPLER.
% QUELLO ORIGINALE E' NELLA DIRECTORY SIMPULSE

%('ENTRA IN SELEZIONE_SKY')
if ~exist('win','var')
    win=1.5;
end
%comm
%  fid=fopen(listin,'r');
% if fid < 0
%     disp([listin ' file could not be opened'])
%     return;
% end
%  nfiles=0;
% 
% while (feof(fid) ~= 1)
%     nfiles=nfiles+1;
%     file{nfiles}=fscanf(fid,'%s',1);
%     str=sprintf('  %s ',file{nfiles});
%     disp(str);
% end
%fine comm
NN=0;

% e=[];
% t=e;
% a=e;
% d=e;

itot=0;
% 
 %for i = 1 : nfiles-1 %coomm
   %fil=file{i}; % commentare 
    %loafil=['load ''' fil ''''];
     %eval(loafil)
     %e=pdata_red{1};    %energia
     len=length(e);
     %disp(sprintf('file %s %d max %e',fil,len,max(e)));
   % N(i)=len;%
    %NN=NN+N(i);%
     %a=pdata{2};
     %d=pdata{3};
     
     %a=pdata_red{2};
     %d=pdata_red{3};
    ii=find(a > ra0(1) & a < ra0(2));
    iii=find(d(ii) > dec0(1) & d(ii) < dec0(2));
    iii=ii(iii);
    itot1=itot+length(iii);
%     for j = 1:22
%         pp=pdata{j};
%         A(itot+1:itot1,j)=pp(iii);
%     end
    itot=itot1;
 %end

NN;

% mua=mean(A(:,2));
% sda=std(A(:,2));
% mud=mean(A(:,3));
% sdd=std(A(:,3));
% 
% in1(1,:)=[mua mud];
% in2(:,1)=A(:,2);
% in2(:,2)=A(:,3);
% ang=astro_angle_m(in1,in2);
% %figure,hist(ang,100)
% 
% ii=find(ang < win*mean(ang));


% MODIFICA MIA
%  ra=in2(ii,1);
%  dec=in2(ii,2);
% figure;plot (ra,dec,'.');

% a=a(iii);
% d=d(iii);
% whos a ;
% whos d;
% pos=[a d];
% whos pos;
% 
% figure, hist3(pos,[50 50])
%figure, surf(pos,[50 50])


%figure; histmat (a,d);
 %figure;(a,d);
% title('AR +/- 7.5 gradi');
% figure;hist(d,50);
% tilte('DEC +/- 7.5 gradi')
