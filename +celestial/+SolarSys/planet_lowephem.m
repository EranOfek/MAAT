function [Pos,Dist,Ang,Mag]=planet_lowephem(JD,Planet,Observer,CooType,Equinox)
% Low-accuracy ephemeris of the planets (~1 arcmin)
% Package: celestial.SolarSys
% Description: Low accuracy ephemeris of the planets.
%                              Accurate to about ~1'.
% Input  : - If single column vector then Julian day.
%            If four column vector then [D M Y FracDay].
%            If six column vector then [D M Y Hour Min Sec].
%          - Planet name:
%            {'Mercury' | 'Venus' | 'Earth' | 'Mars' | 'Jupiter' |
%             'Saturn'  | 'Uranus' | 'Neptune' | 'Sun'}
%          - Observer position within the solar system:
%            {'Mercury' | 'Venus' | 'Earth' | 'Mars' | 'Jupiter' |
%             'Saturn'  | 'Uranus' | 'Neptune' | 'Sun'},
%            default is 'Earth' (geocentric).
%          - Coordinate type:
%            'RectEcl'    - Rectangular ecliptic [au, au, au]
%            'RectEq'     - Rectangular equatorial [au, au, au]
%            'SphericEcl' - Spherical ecliptic [rad, rad, au]
%            'SphericEq'  - Spherical equatorial [rad, rad, au] - default.
%          - Equinox :
%            'J2000' (not supported in this version)
%            'date'   - default
% Output : - Matrix of coordinates [X, Y, Z] or [RA, Dec, RadVec]
%            row per each date.
%          - Matrix of distances:
%            [Delta(observer-object) RadiusVector(Sun-Object), (Sun-Observer)]
%          - Matrix of angles:
%            [Sun-Observer-Planet,  Sun-Planet-Observer] in radians.
%          - [Magnitude, IlluminatedFraction]
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Pos,Dist,Ang,Mag]=celestial.SolarSys.planet_lowephem(2451545+(0:1:100)','Mars','Earth','SphericEq','date');
% Reliable: 2
%--------------------------------------------------------------------------
RAD  = 180./pi;

if (nargin==2),
   Observer = 'Earth';
   CooType  = 'SphericEq';
   Equinox  = 'date';
elseif (nargin==3),
   CooType  = 'SphericEq';
   Equinox  = 'date';
elseif (nargin==4),
   Equinox  = 'date';
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end


if (size(JD,2)==1),
   % do nothing - allready in JD
elseif (size(JD,2)==4),
   JD = julday(JD);
elseif (size(JD,2)==6),
   FracDay = celestial.coo.convertdms(JD(:,4:6),'H','f');
   JD      = celestial.time.julday([JD(:,1:3), FracDay]);
else
   error('Illegal number of columns in JD input');
end

Nt = length(JD);
%--- Observer coordinates ---
switch Observer
 case 'Sun'
    ObserverCoo = zeros(Nt,3);
 otherwise
    % rest of planets
    [X, Y, Z]   = celestial.SolarSys.ple_planet(JD,Observer,'XYZ');
    ObserverCoo = [X, Y, Z];
end

%--- Planet coordinates ---
switch Planet
 case 'Sun'
    PlanetCoo   = zeros(Nt,3);
 otherwise
    % rest of planets
    [X, Y, Z]   = celestial.SolarSys.ple_planet(JD,Planet,'XYZ');
    PlanetCoo   = [X, Y, Z];
end


PosRectEcl = PlanetCoo - ObserverCoo;
% [Delta, R, Sun-Observer]
Dist = [sqrt(sum(PosRectEcl.^2, 2)),...
        sqrt(sum(PlanetCoo.^2, 2)),...
        sqrt(sum(ObserverCoo.^2, 2))];
% [Sun-Observer-Planet,  Sun-Planet-Observer]
Ang  = [acos((Dist(:,1).^2 + Dist(:,3).^2 - Dist(:,2).^2)./(2.*Dist(:,1).*Dist(:,3))),...
        acos((Dist(:,1).^2 + Dist(:,2).^2 - Dist(:,3).^2)./(2.*Dist(:,1).*Dist(:,2)))];
% [Magnitude, Illuminated fraction]
Mag = celestial.SolarSys.planets_magnitude(Planet, Dist(:,2), Dist(:,1), Dist(:,3),'V',0,0);
K   = 0.5.*(1 + cos(Ang(:,2)));
Mag = [Mag, K];

switch CooType
 case 'RectEcl'
    Pos = PosRectEcl;
    %--- Equinox ---
    switch Equinox
     case 'date'
        %--- allready relative to date ---
     case 'J2000'
        error('J2000 coordinates are not available in RectEcl mode');
     otherwise
        error('Unknown Equinox option');
    end
 case 'SphericEcl'
    Pos = [atan2(PosRectEcl(:,2),PosRectEcl(:,1)),...
           atan(PosRectEcl(:,3)./sqrt( PosRectEcl(:,2).^2 + PosRectEcl(:,1).^2 ) ),...
           Dist(:,1)];
    %--- Equinox ---
    switch Equinox
     case 'date'
        %--- allready relative to date ---
     case 'J2000'
        error('J2000 coordinates are not available in RectEcl mode');
     otherwise
        error('Unknown Equinox option');
    end
 case {'RectEq','SphericEq'}
    switch Equinox
     case 'J2000'
        error('J2000 not supported in this version');
    end

    PosRectEq = zeros(Nt,3);
    for I=1:1:Nt,
       % ecliptic with true ecliptic and equinox of date
       % to equatorial with true equinox of date
       RotMat         = celestial.coo.rotm_coo('E',JD(I));
       PosRectEq(I,:) = [RotMat*PosRectEcl(I,:).'].';
    end

    switch CooType
     case 'RectEq'
        Pos = PosRectEq;
     case 'SphericEq'
        Pos = [atan2(PosRectEq(:,2),PosRectEq(:,1)),...
               atan(PosRectEq(:,3)./sqrt( PosRectEq(:,2).^2 + PosRectEq(:,1).^2 ) ),...
               Dist(:,1)];
     otherwise
        error('Unknown CooType option');
    end
 otherwise
    error('Unknown CooType option');
end
