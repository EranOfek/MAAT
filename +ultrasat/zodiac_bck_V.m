function [ZodiVmag]=zodiac_bck_V(Long,Lat,Date,varargin)
% Get the Zodiac V light surface brightness as a function of coordinates
% Package: @AstSpec
% Description: Calculate the zodiac V vand Vega magnitude for a sky.
% Input  : - Either Ecliptic longitude or Helio-ecliptic longitude
%            (i.e., L - L_sun) in radians
%            (see convertdm.m for additional options).
%            If date is not given, this is assumed to be Helio-ecliptic
%            longitude.
%          - Ecliptic latitude [radians].
%            (see convertdm.m for additional options).
%          - Date [D M Y] or JD at which to calculate the Helio-ecliptic
%            longitude. Defaut is empty. If empty than treat the first
%            input argument as Helio-ecliptic longitude.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'InterpMethod' - Default is 'linear'.
%            'ReplaceNaNVal' - Value to replace NaNs. Default is 21.5.
% Output : - Zodiacal light V-vabd Vega surface magnitude.
% Tested : Matlab R2014a
%     By : Ilan Sagiv                      Sep 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: https://hst-docs.stsci.edu/display/WFC3IHB/9.7+Sky+Background#id-9.7SkyBackground-9.7.19.7.1ZodiacalLight,EarthShine,andLOW-SKY
% Example: [ZodiVmag]=ultrasat.zodiac_bck_V(Long,Lat,Date)
% Reliable: 2
%--------------------------------------------------------------------------

RAD            = 180./pi;
h              = constant.h; %('h','cgs');
c              = constant.c; %('c','cgs');
hc             = h*c*1e8; %[erg]*[Ang]

if (nargin<3)
    Date = [];
end



DefV.InterpMethod         = 'linear';
DefV.ReplaceNaNVal        = 21.5;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (~isempty(Date))
    % assume ecliptic longitude as input
    % Transform Ecliptic coordinates to Helio-ecliptic coordinates
    HelioEcLong = celestial.coo.ecliptic2helioecliptic(Long,Date);
else
    HelioEcLong = Long;
end

% convert to 0-pi range
HelioEcLong(HelioEcLong>pi) = 2.*pi - HelioEcLong(HelioEcLong>pi);


% get sodi spectrum
Spec = ultrasat.zodiac_spectrum;

% approximate zodiacal sky background as a function of
% Heliocentric ecliptic longitude and ecliptic latitude
% (in Vmag / arcsec^2)
% from HST WFC3 Handbook; (table 9.4)
% V-band magnitude of zodi
% adopted from Table 9.4 in:
% http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c09_exposuretime08.html#389841
EcLon=(0:15:180)./RAD;
EcLat=(0:15:90)./RAD;
%lat-->
%v lon
ZodiVmagTable=[...
nan ,nan ,nan ,nan ,22.6,23.0,23.3;
nan ,nan ,nan ,nan ,22.6,23.1,23.3;
nan ,nan ,nan ,22.3,22.7,23.1,23.3;
nan ,nan ,22.1,22.5,22.9,23.1,23.3;
21.3,21.9,22.4,22.7,23.0,23.2,23.3;
21.7,22.2,22.6,23.9,23.1,23.2,23.3;
22.0,22.3,22.7,23.0,23.2,23.3,23.3;
22.2,22.5,22.9,23.1,23.3,23.3,23.3;
22.4,22.6,22.9,23.2,23.3,23.3,23.3;
22.4,22.6,22.9,23.2,23.3,23.3,23.3;
22.4,22.6,22.9,23.1,23.3,23.3,23.3;
22.3,22.5,22.8,23.0,23.2,23.4,23.3;
22.1,22.4,22.7,23.0,23.2,23.4,23.3];

% Interpolate zodi Vmag to requested positions
ZodiVmag = interp2(EcLat,EcLon,ZodiVmagTable,abs(Lat),HelioEcLong,InPar.InterpMethod);
ZodiVmag(isnan(ZodiVmag)) = InPar.ReplaceNaNVal;
