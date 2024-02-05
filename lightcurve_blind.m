function [lc]=lightcurve_blind(t,per,np,x)

% p_lightcurve  light curve for pulses
%
%     lc=p_lightcurve(t,per,np)
%
%   t     pulses
%   per   period
%   np    number of bins in the period

tt=mod(t,per);
% dx=per/np;
% x=dx/2:dx:per;
lc=hist(tt,x);
lc(np+1)=lc(1);

% lc=gd(lc);
% lc=edit_gd(lc,'dx',dx);