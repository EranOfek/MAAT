function M=planets_magnitude(PlanetName, R, Delta, Re, Band, DU, B)
% Calculate the planets apparent magnitude
% Package: celestial.SolarSys
% Description: Calculate the planets apparent magnitude.
% Input  : - Planet name:
%            'Mercury' | 'Venus' | 'Earth' | 'Mars' | 'Jupiter' |
%            'Saturn' | 'Uranus' | 'Neptune' | 'Pluto'
%          - Sun-planet distance in au.
%          - Observer-planet distance in au.
%          - Sun observer distance in au.
%          - Band (filter) type:
%            'V' : visual (default)
%            'U' : UV
%            'B' : blue
%          - The difference between the Saturnicentric longitudes of the Sun
%            and the Earth (radians),
%            relavant only for Saturn.
%          - Saturnocentric sub-Earth latitude (radians),
%            relavant only for Saturn.
% Output : - The planet magnitude
% Reference: A&A supp.
% Tested : matlab 5.3
%     By : Eran O. Ofek                    May 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: M=celestial.SolarSys.planets_magnitude('Mars',1.5,0.6,1,'V');
% Reliable: 2
%--------------------------------------------------------------------------
RADIAN = 180./pi;

if (nargin==4),
   Band = 'V';
   DU = NaN;
   B = NaN;
elseif (nargin==5),
   DU = NaN;
   B = NaN;
elseif (nargin==6),
   B = NaN;
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end 
DU = DU.*RADIAN;
B  = B.*RADIAN;


% the phase angle
I = acos((R.^2 + Delta.^2 - Re.^2)./(2.*R.*Delta));
I = I.*RADIAN;  % convert to deg.

switch Band
 case 'V'
    BandOffset = 0;
 case 'U'
    switch lower(PlanetName)
     case 'mercury'
        BandOffset = 0.93 + 0.41;
     case 'venus'
        BandOffset = 0.82 + 0.50;
     case 'earth'
        BandOffset = NaN;
     case 'mars'
        BandOffset = 1.36 + 0.58;
     case 'jupiter'
        BandOffset = 0.83 + 0.48;
     case 'saturn'
        BandOffset = 1.04 + 0.58;
     case 'uranus'
        BandOffset = 0.56 + 0.28;
     case 'neptune'
        BandOffset = 0.41 + 0.21;
     case 'pluto'
        BandOffset = 0.80 + 0.31;
     otherwise
        error('Unknown planet');
    end
 case 'B'
    switch PlanetName
     case 'Mercury'
        BandOffset = 0.93;
     case 'Venus'
        BandOffset = 0.82;
     case 'Earth'
        BandOffset = NaN;
     case 'Mars'
        BandOffset = 1.36;
     case 'Jupiter'
        BandOffset = 0.83;
     case 'Saturn'
        BandOffset = 1.04;
     case 'Uranus'
        BandOffset = 0.56;
     case 'Neptune'
        BandOffset = 0.41;
     case 'Pluto'
        BandOffset = 0.80;
     otherwise
        error('Unknown planet');
    end
 otherwise
    error('Unknown band');
end



switch lower(PlanetName)
 case 'mercury'
    M = -0.42 + 5.*log10(R.*Delta) + 0.0380.*I - 0.000273.*I.^2 + 0.000002.*I.^3;
 case 'venus'
    M = -4.40 + 5.*log10(R.*Delta) + 0.0009.*I + 0.000239.*I.^2 - 0.00000065.*I.^3;
 case 'earth'
    M = -3.86 + 5.*log10(R.*Delta);
 case 'mars'
    M = -1.52 + 5.*log10(R.*Delta) + 0.016.*I;
 case 'jupiter'
    M = -9.40 + 5.*log10(R.*Delta) + 0.005.*I;
 case 'saturn'
 M = -8.88 + 5.*log10(R.*Delta) + 0.044.*abs(DU) - 2.60.*sin(abs(B)) + 1.25.*sin(abs(B)).^2; 
 case 'uranus'
    M = -7.19 + 5.*log10(R.*Delta);
 case 'neptune'
    M = -6.87 + 5.*log10(R.*Delta);
 case 'pluto'
    M = -1.00 + 5.*log10(R.*Delta);
 otherwise
    error('Unknown planet');
end


M = M + BandOffset;
 
