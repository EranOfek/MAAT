function [Alt,AM]=ha2alt(HA,Dec,Lat)
% Hour angle to altitude and airmass
% Package: celestial.coo
% Description: Given Hour Angle as measured from the meridian, the source
%              declination and the observer Geodetic latitude, calculate
%              the source altitude above the horizon and its airmass.
% Input  : - Hour Angle [radians].
%          - Declination [radians].
%          - Latitude [radians].
% Output : - Altitude [radians].
%          - Airmass.
% See also: horiz_coo.m, ha2az.m
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Aug 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [A,AM]=celestial.coo.ha2alt(1,1,0.5)
% Reliable: 1
%--------------------------------------------------------------------------

Alt = asin(sin(Dec).*sin(Lat) + cos(Dec).*cos(Lat).*cos(HA));
if (nargout>1),
    AM  = celestial.coo.hardie(pi./2-Alt);
end
