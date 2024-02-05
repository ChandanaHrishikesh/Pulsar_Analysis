function [ra dec]=hour2deg_coord(ra,dec) % from hours to deg
ra=hour2deg(ra);


if dec(1)<0
   dec=dec(1)-dec(2)/60-dec(3)/3600;
else
    dec=dec(1)+dec(2)/60+dec(3)/3600;
end
