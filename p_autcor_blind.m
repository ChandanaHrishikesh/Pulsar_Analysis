function [autc sp_re sp_im]=p_autcor_blind(t,T,dt)
% p_autcor  impulsive autocorrelation
%
%      [autc sp per]=p_autcor(t,T,dt)
%
%   autc    impulsive autocorrelation
%   sp      spectrum
%
%   t       pulses
%   T       max delay
%   dt      bin width

T= T*86400;

n=length(t);            % nume di fotoni
N=round(T/dt)           % numero di bin all'interno di un Tp

% modifica**************************
% NFFT = 2^nextpow2(N);
% dt=T/NFFT;
% N=NFFT;
%**********************************


autc=zeros(1,N);        % autocorrelazione

for i = 1:n-1
  
    k=ceil((t(i+1:n)-t(i))/dt+0.5);
    ii=find(k<=N);
    k=k(ii);
    autc(k)=autc(k)+1;
end

autc(1)=autc(1)*2+n;
autc1=zeros(1,2*N);
autc1(1:N)=autc;
autc1(2*N:-1:N+2)=autc(2:N);
autc1(N+1)=autc1(N);


sp=fft(autc1)*dt;

% MODIFICA DEL 20 DICEMBRE 2012********************************************
% FACCIO LO SPETTRO SU META' INTERVALLO

sp=sp(2:round(length(sp)/2));

%*************************************************************************






autc=gd(autc);
autc=edit_gd(autc,'dx',dt);

sp_re=gd(real(sp));
sp_im=gd(imag(sp));

sp_re=edit_gd(sp_re,'dx',1/(2*N*dt));
sp_im= edit_gd(sp_im,'dx',1/(2*N*dt));