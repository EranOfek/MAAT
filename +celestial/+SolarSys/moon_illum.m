function [Illum,Ph]=moon_illum(Date)
% Low accuracy Moon illuminated fraction
% Package: celestial.SolarSys
% Description: Low accuracy Moon illuminated fraction
% Input  : - JD or date (see julday.m for options).
% Output : - Illuminated fraction of the moon.
%          - Moon Phase angle [radians].
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Aug 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Illum,Ph]=celestial.SolarSys.moon_illum(2451545)
% Reliable: 1
%------------------------------------------------------------------------------
RAD = 180./pi;

if (size(Date,2)==1),
   % already JD
   JD = Date;
else
   JD = celestial.time.julday(Date);
end

T   = (JD - 2451545.0)./36525;

D  = 297.8502042 + 445267.1115168.*T - 0.0016300.*T.*T + T.*T.*T./545868 - T.^4./113065000;
Mt = 134.9634114 + 477198.8676313.*T + 0.0089970.*T.*T + T.*T.*T./69699 - T.^4./14712000;
M  = 357.5291092 + 35999.0502909.*T - 0.0001536.*T.*T + T.*T.*T./24490000;

D  = D./RAD;
Mt = Mt./RAD;
M  = M./RAD;

Ph  = 180 - D.*RAD - 6.289.*sin(Mt) + ...
                     2.100.*sin(M) - ...
                     1.274.*sin(2.*D - Mt) - ...
                     0.658.*sin(2.*D) - ...
                     0.214.*sin(2.*Mt) - ...
                     0.110.*sin(D);

Ph = Ph./RAD;
Ph = (Ph./(2.*pi) - floor(Ph./(2.*pi))).*2.*pi;
Illum = 0.5.*(1+cos(Ph));

if (Ph>pi),
   Illum = -Illum;
end
