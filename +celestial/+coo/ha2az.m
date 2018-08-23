function [Az,Alt,AM]=ha2az(HA,Dec,Lat)
% Convert hour angle and declination to azimuth, altitude and airmass
% Package: celestial.coo
% Description: Given Hour Angle as measured from the meridian, the source
%              declination and the observer Geodetic latitude, calculate
%              the horizonal source azimuth
% Input  : - Hour Angle [radians].
%          - Declination [radians].
%          - Latitude [radians].
% Output : - Azimuth [radians].
%          - Altitude [radians].
%          - Airmass.
% See also: horiz_coo.m, ha2alt.m
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Aug 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Az,Alt,AM]=celestial.coo.ha2az(1,1,1)
% Reliable: 1
%--------------------------------------------------------------------------


SinAlt = sin(Dec).*sin(Lat) + cos(Dec).*cos(HA).*cos(Lat);
CosAlt = sqrt(1-SinAlt.*SinAlt);

SinAz  = (-cos(Dec).*sin(HA))./CosAlt;
CosAz  = (sin(Dec).*cos(Lat) - cos(Dec).*cos(HA).*sin(Lat))./CosAlt;

Az     = atan2(SinAz, CosAz);
if (nargin>1),
    Alt = asin(sin(Dec).*sin(Lat) + cos(Dec).*cos(Lat).*cos(HA));
    if (nargin>2),
        AM  = celestial.coo.hardie(pi./2-Alt);
    end
end
