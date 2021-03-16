function [Az,Alt]=hadec2azalt(HA,Dec,Lat,CooUnits)
% Convert HA/Dec to Az/Alt
% Package: +celestial.coo
% Input  : - HA [rad]
%          - Dec [rad]
%          - Lat [rad]
%          - Input/Output coordinates units {'rad' | 'deg'}. Default is 'rad'.
% Output : - Az [same as input]
%          - Alt [same as input]
%      By : Eran O. Ofek                Aug 2020
% Example: [Az,Alt]=celestial.coo.hadec2azalt(HA,Dec,Lat)

if nargin==4
    Convert = convert.angular(CooUnits,'rad');
    HA  = HA.*Convert;
    Dec = Dec.*Convert;
    Lat = Lat.*Convert;
end


SinAlt = sin(Dec).*sin(Lat) + cos(Dec).*cos(HA).*cos(Lat);
CosAlt = sqrt(1-SinAlt.*SinAlt);
SinAz  = (-cos(Dec).*sin(HA))./CosAlt;
CosAz  = (sin(Dec).*cos(Lat) - cos(Dec).*cos(HA).*sin(Lat))./CosAlt;
Az     = atan2(SinAz, CosAz);
Alt    = asin(SinAlt);
I      = Az<0;
Az(I)  = 2.*pi+Az(I);

if nargin==4
    Convert = convert.angular('rad',CooUnits);
    Az  = Az.*Convert;
    Alt = Alt.*Convert;
end
