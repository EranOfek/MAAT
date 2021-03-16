function [HA,Dec]=azalt2hadec(Az,Alt,Lat,CooUnits)
% Convert Az/Alt to HA/Dec
% Package: +celestial.coo
% Input  : - Az [rad]
%          - Alt [rad]
%          - Lat [rad]
%          - CooUnits {'rad' | 'deg'}. Default is 'rad'.
% Output : - HA [same as input]
%          - Dec [same as input]
%      By : Eran O. Ofek                Aug 2020
% Example: Az=1; Alt=1; Lat=[30;31]; [HA,Dec]=celestial.coo.azalt2hadec(Az,Alt,Lat,'deg')

if nargin==4
    Convert = convert.angular(CooUnits,'rad');
    Az  = Az.*Convert;
    Alt = Alt.*Convert;
    Lat = Lat.*Convert;
end

SinDec = sin(Alt).*sin(Lat) + cos(Alt).*cos(Az).*cos(Lat);
CosDec = sqrt(1 - SinDec.*SinDec);

SinHA  = (-cos(Alt).*sin(Az))./CosDec;
CosHA  = (sin(Alt).*cos(Lat) - cos(Alt).*cos(Az).*sin(Lat))./CosDec;
HA     = atan2(SinHA, CosHA);
Dec    = asin(SinDec);


if nargin==4
    Convert = convert.angular('rad',CooUnits);
    HA  = HA.*Convert;
    Dec = Dec.*Convert;
end
