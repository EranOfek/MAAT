function [Alpha0,Delta0,W]=planets_rotation(PlanetName,JD,JupiterSystem)
% Planet north pole, rotation rate and the primery meridian
% Package: celestial.SolarSys
% Description: Return Planet north pole, rotation rate and the primery
%              meridian.
% Input  : - PlanetName {'Sun' | 'Mercury' | 'Venus' | 'Earth' |
%                        'Mars' | 'Jupiter' | 'Saturn' | 'Uranus' |
%                        'Neptune' | 'Pluto'}
%          - Vector of Julian days.
%          - Jupiter rotation system {'I' | 'II' | 'III'},
%            default is 'I'.
% Output : - Alpha0 : J2000.0 Right Ascension of the north pole [radians].
%          - Delta0 : J2000.0 Declination of the north pole [radians].
%          - W : location of the prime meridian measured along the
%                planet's equator in an easterly direction with
%                respect to the planet's north pole from the node
%                (located at RA 90deg+Alpha0) of the planet's
%                equator on the standard equator [radians].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Dec 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Alpha0,Delta0,W]=celestial.SolarSys.planets_rotation('Jupiter',2451545,'II')
% Reliable: 2
%--------------------------------------------------------------------------

RAD          = 180./pi;
JD2000        = 2451545.0;
DaysFromJ2000 = JD - JD2000;
T             = DaysFromJ2000./36525;

if (nargin==2),
   JupiterSystem = 'I';
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end


switch PlanetName
 case 'Sun'
    Alpha0  = 286.13;
    Delta0  =  63.87;
    W       = 84.10 + 14.1844000.*DaysFromJ2000;
 case 'Mercury'
    Alpha0  = 281.01 - 0.003.*T;
    Delta0  = 61.45  - 0.005.*T;
    W       = 329.71 + 6.1385025.*DaysFromJ2000;
 case 'Venus'
    Alpha0  = 272.72;
    Delta0  =  67.15;
    W       = 160.26 - 1.4813596.*DaysFromJ2000;
 case 'Earth'
    Alpha0  = 0.00  - 0.641.*T;
    Delta0  = 90.00 - 0.557.*T;
    W       = 190.16 + 360.9856235.*DaysFromJ2000;
 case 'Mars'
    Alpha0  = 317.681 - 0.108.*T;
    Delta0  =  52.886 - 0.061.*T;
    W       = 176.868 + 350.8919830.*DaysFromJ2000;
 case 'Jupiter'
    Alpha0  = 268.05 - 0.009.*T;
    Delta0  =  64.49 + 0.003.*T;

    switch JupiterSystem
     case 'III'
        % rotation of Jupiter magnetic field
        W      = 284.95 + 870.5360000.*DaysFromJ2000;
     case 'II'
        % mean atmospheric rotation north of the south
        % componenet of the north equatorial belt,
        % and south of the north component of the south
        % equatorial belt.
        W      = 43.3  + 870.270.*DaysFromJ2000;
     case 'I'
        % mean atmospheric equatorial rotation
        W      = 67.1  + 877.900.*DaysFromJ2000;
     otherwise
        error('Unknown JupiterSystem option');
    end
 case 'Saturn'
    Alpha0  = 40.58 - 0.036.*T;
    Delta0  = 83.54 - 0.004.*T;
    % magnetic field rotation
    W       = 38.90 + 810.7939024.*DaysFromJ2000;
 case 'Uranus'
    Alpha0  = 257.43;
    Delta0  = -15.10;
    % magnetic field rotation
    W       = 203.81 - 501.1600928.*DaysFromJ2000;
 case 'Neptune'
    N       = (359.28 + 54.308.*T)./RAD;
    Alpha0  = 299.36 + 0.70.*sin(N);
    Delta0  =  43.46 - 0.51.*cos(N);
    % magnetic field rotation
    W       = 253.18 + 536.3128492.*DaysFromJ2000 - 0.48.*sin(N);
 case 'Pluto'
    Alpha0  = 313.02;
    Delta0  =   9.09;
    W       = 236.77 - 56.3623195.*DaysFromJ2000;
 otherwise
    error('Unknwon PlanetName option');
end



Alpha0 = Alpha0./RAD;
Delta0 = Delta0./RAD;
W      = W./RAD;
