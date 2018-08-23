function [Phys, Ang]=planet_radius(PlanetName,Distance,PlLat)
% Planet radius and flattening factor and angular diameter.
% Package: celestial.SolarSys
% Description: Get planetary radius and flattening factor, and calculate
%              planet angular diameter.
% Input  : - Planet name:
%            'Mercury' | 'Venus' | 'Earth' | 'Mars' | 'Jupiter' |
%            'Saturn' | 'Uranus' | 'Neptune' | 'Pluto' | 'Sun' | 'Moon'
%          - (optional) Observer-planet distance [au].
%            If the distance is given, the angular diameter is also
%            calculated.
%          - Planetocentric latitude of the observer, Latitude
%            of subobserver [radians]. (default is 0).
% Output : - [Planetary equatorial radius (km),
%             Planet polar radius (km),
%             Flattening factor]
%          - [Angular equatorial radius, (radians)
%             Angular polar radius, (radians)],
%             The angular equatorial radius is given only if the distance
%             is given, and the ang. polar radius is given only if the
%             latitude sub-observer is given.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    May 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%----------------------------------------------------------------------------
AU = 149597870.66; %km

switch PlanetName
 case 'Mercury'
    A = 2439.7;
    F = 0;
 case 'Venus'
    A = 6051.9;
    F = 0;
 case 'Earth'
    A = 6378.140;
    F = 1./298.257
 case 'Moon'
    A = 1738;
    F = 0;
 case 'Mars'
    A = 3397;
    F = 0.00647630;
 case 'Jupiter'
    A = 71492;
    F = 0.0648744;
 case 'Saturn'
    A = 60268;
    F = 0.0979624;
 case 'Uranus'
    A = 25559;
    F = 0.0229273;
 case 'Neptune'
    A = 24764;
    F = 0.0171;
 case 'Pluto'
    A = 1151;
    F = 0;
 case 'Sun'
    A = 696000;
    F = 0;
 otherwise
    error('Unknown object name');
end

B = A - F.*A;

if (nargin==1),
   Phys = [A, B, F];
   Ang = NaN;   
elseif (nargin==2 | nargin==3),
   Aa    = atan(A./(AU.*Distance));
   if (nargin==2),
      Gamma = zeros(size(Distance));
   end
   Gamma = pi./2 - abs(PlLat);
   SinG2 = sin(Gamma).^2;
   F2    = F.^2;
   Ba = Aa.*(1 - 0.5.*(2.*F.*SinG2 - F2.*SinG2) + ((2.*F.*SinG2 - F2.*SinG2).^2)./8);
   Ang  = [Aa, Ba];   
   Phys = [A, B, F];
else
   error('Illegal number of input arguments');
end
