function [Tot_Illum, Sun_Illum, Moon_Illum]=skylight(JD,GeodPos)
% Calculate the total sky illumination in Lux on horizontal surface
% Package: celestial.SolarSys
% Description: Calculate the total sky illumination due to the Sun, Moon,
%              stars and air-glow, in Lux on horizontal surface as a 
%              function of time.
% Input  : - vector of JD.
%          - Geodetic position [East_Long, North_Lat] in radians.
% Output : - Total illumination in Lux on horiz. surface.
%          - Sun+sky illumination in Lux on horiz. surface.
%          - Moon illumination in Lux on horiz. surface.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Tot_Illum, Sun_Illum, Moon_Illum]=skylight(2451545,[1 1])
% Reliable: 2
%------------------------------------------------------------------------------

[Sun_RA, Sun_Dec]            = celestial.SolarSys.suncoo(JD);
[Moon_RA, Moon_Dec, Moon_HP] = celestial.SolarSys.mooncool(JD,GeodPos);

Moon_Elon = acos(sin(Sun_Dec).*sin(Moon_Dec) + cos(Sun_Dec).*cos(Moon_Dec).*cos(Sun_RA-Moon_RA));

Sun_Horiz  = celestial.coo.horiz_coo([Sun_RA,Sun_Dec],JD,GeodPos,'h');
Moon_Horiz = celestial.coo.horiz_coo([Moon_RA,Moon_Dec],JD,GeodPos,'h');
% calculate sunlight + starlight + airglow
Sun_Illum  = celestial.SolarSys.sunlight(Sun_Horiz(:,2));

% calculate moonlight
Moon_Illum = celestial.SolarSys.moonlight(Moon_Horiz(:,2),Moon_HP,Moon_Elon);

Tot_Illum = Sun_Illum + Moon_Illum;

