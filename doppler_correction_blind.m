function [t_corr]=doppler_correction_blind (t_fotoni_mjd,n);




format long;
%KMAU= 1.49597870691000015e+08;  % conversione Km-AU
c_light= 299792.458;            % velocita' della luce in Km/s

rif_gps=mjd2gps(51910);               % mjd  di rifermento convertito in gps




[p_dec_terra p_dec_sc]=pos_detector_blind(t_fotoni_mjd);

%l=length(t);


x_terra=p_dec_terra(1,:);
y_terra=p_dec_terra(2,:);
z_terra=p_dec_terra(3,:);


x_sc=p_dec_sc(1,:);
y_sc=p_dec_sc(2,:);
z_sc=p_dec_sc(3,:);


x_sc_ssb= x_sc + x_terra;     %posizione x SC-SSB
y_sc_ssb= y_sc + y_terra;
z_sc_ssb= z_sc + z_terra;
 

        
posx=x_sc_ssb*n(1);
posy=y_sc_ssb*n(2);
posz=z_sc_ssb*n(3);

       
somma=(posx+posy+posz)/c_light;

t_fotoni_gps=mjd2gps(t_fotoni_mjd);  


t_corr=somma+t_fotoni_gps;

[shapiro_delay]= correzione_shapiro (n, x_sc_ssb, y_sc_ssb, z_sc_ssb);

t_corr=t_corr- shapiro_delay;


%  for jj=1:length(t_fotoni_gps)
%      fprintf(tempi,'%15.10f\n',t_corr(jj));
%  end
    