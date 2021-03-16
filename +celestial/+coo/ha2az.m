function [Az,Alt,AM]=ha2az(HA,Dec,Lat,Units)
% Convert hour angle and declination to azimuth, altitude and airmass
% Package: celestial.coo
% Description: Given Hour Angle as measured from the meridian, the source
%              declination and the observer Geodetic latitude, calculate
%              the horizonal source azimuth
% Input  : - Hour Angle [radians].
%          - Declination [radians].
%          - Latitude [radians].
%          - Input/output units. 'rad' | 'deg'. Default is 'rad'.
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


if nargin<4
    Units = 'rad';
end

switch lower(Units)
    case 'rad'
        % do nothing
    otherwise
        Convert = convert.angular(Units,'rad');
        HA      = HA.*Convert;
        Dec     = Dec.*Convert;
        Lat     = Lat.*Convert;
end
          

SinAlt = sin(Dec).*sin(Lat) + cos(Dec).*cos(HA).*cos(Lat);
CosAlt = sqrt(1-SinAlt.*SinAlt);

SinAz  = (-cos(Dec).*sin(HA))./CosAlt;
CosAz  = (sin(Dec).*cos(Lat) - cos(Dec).*cos(HA).*sin(Lat))./CosAlt;

Az      = atan2(SinAz, CosAz);
Convert = convert.angular('rad',Units);

if (nargin>1)
    Alt = asin(sin(Dec).*sin(Lat) + cos(Dec).*cos(Lat).*cos(HA));
    if (nargin>2)
        AM  = celestial.coo.hardie(pi./2-Alt);
    end
    Alt = Alt.*Convert;
end

Az = Az.*Convert;
