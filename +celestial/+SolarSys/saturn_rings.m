function [SLO,SLS,DU,P,Major,Minor]=saturn_rings(JD)
% Calculate the orientation angles for Saturn's rings
% Package: celestial.SolarSys
% Description: Calculate the orientation angles for Saturn's rings.
% Input  : - Time of observation, JD, TT time scale.
% Output : - The saturnicentric latitude of the observer referred to the
%            plane of the ring, positive towereds the north, [radians].
%          - Saturnicentric latitude of the Sun referred to the
%            plane of the ring; positive towerds the north, [radians].
%          - Difference between saturnicentric longitude of sun
%            and observer, [radians].
%          - The geocentric position angle of the northern semi minor axis
%            of the apparent ellipse of the ring, measured from the north
%            towerds the east, [radians].
%          - Major axis of ring [arcsec].
%            Columns 1 to 5 are for:
%            Outer edge of outer ring
%            Inner edge of outer ring
%            Outer edge of inner ring
%            Inner edge of inner ring
%            Inner edge of dusky ring
%          - Minor axis of ring [arcsec].
%            Columns 1 to 5 are for:
%            Outer edge of outer ring
%            Inner edge of outer ring
%            Outer edge of inner ring
%            Inner edge of inner ring
%            Inner edge of dusky ring
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    May 2001
% Web example: http://astroclub.tau.ac.il/ephem/Saturn/index.php
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [SLO,SLS,DU,P,Major,Minor]=celestial.SolarSys.saturn_rings(2451545);
% Reliable: 2
%---------------------------------------------------------------------------
RAD = 180./pi;
if (nargin==1),
   [Coo,Dist]=celestial.SolarSys.planet_ephem(JD,'Saturn','Earth','ecliptic');
   R     = Dist(:,2);
   Delta = Dist(:,1);
   L     = Coo(:,1);
   B     = Coo(:,2);
   [Coo,Vel] = celestial.SolarSys.calc_vsop87(JD,'Saturn','d','d');
   Ls    = Coo(1,:).';
   Bs    = Coo(2,:).';
else
   error('Illegal number of input arguments');
end

%L*RAD
%B*RAD
%R*RAD
%Delta*RAD

T = (JD - 2451545)./36525;

% Inclination of the plane of rings
% and the longitude of the asc. node
% referred to the ecliptic and mean equinox of B1950.
%I     = 28.0817;  % deg.
%Omega = 168.8112; % deg.

% convert to equinox of date
I     = 28.075216 - 0.012998.*T + 0.000004.*T.^2;
Omega = 169.508470 + 1.394681.*T + 0.000412.*T.^2;

% convert to radians
I     = I./RAD;
Omega = Omega./RAD;

% the saturnicentric latitude of the observer referred to the plane
% of the ring, positive towereds the north.
SLO    =  asin(sin(I).*cos(B).*sin(L-Omega) - cos(I).*sin(B));

% ring factors:
% Outer edge of outer ring
% Inner edge of outer ring
% Outer edge of inner ring
% Inner edge of inner ring
% Inner edge of dusky ring
Factor = [1, 0.8801, 0.8599, 0.6650, 0.5486];
% ring major axis [arcsec.]
Major     = 375.35./Delta;
% ring minor axis [arcsec.]
Minor     = Major.*sin(abs(SLO));

Major     = Major*Factor;
Minor     = Minor*Factor;

% longitude of the ascending node of Saturn's orbit
N = 113.6655 + 0.8771.*T;
N = N./RAD;

% correct Ls and Bs for the Sun's abberation
Bs = Bs - 0.000764.*cos(Ls-N)./(RAD.*R);
Ls = Ls - 0.01759./(RAD.*R);

% Saturnicentric latitude of the Sun referred to the
% plane of the ring; positive towerds the north.
SLS = asin(sin(I).*cos(Bs).*sin(Ls-Omega) - cos(I).*sin(Bs));

% diff. between saturnicentric long of sun and observer
U1 = atan((sin(I).*sin(Bs) + cos(I).*cos(Bs).*sin(Ls-Omega))./(cos(Bs).*cos(Ls-Omega)));
U2 = atan((sin(I).*sin(B) + cos(I).*cos(B).*sin(L-Omega))./(cos(B).*cos(L-Omega)));
DU = abs(U1-U2);

L0 = Omega - pi./2;
L  = L + 0.005693.*cos(L0-L)./(RAD.*cos(B));
B0 = pi./2 - I;
B  = B + 0.005693.*sin(L0-L).*sin(B)./RAD;

[Nut] = celestial.coo.nutation(JD,'f');
L  = L  + Nut(:,1);
L0 = L0 + Nut(:,1);

RA   = zeros(size(L));
Dec  = zeros(size(B));
RA0  = zeros(size(L0));
Dec0 = zeros(size(B0));
for J=1:1:length(JD),
   RotM    = celestial.coo.rotm_coo('Ed',JD(J));
   X       = celestial.coo.cosined([L,B]);
   X0      = celestial.coo.cosined([L0,B0]);
   C       = celestial.coo.cosined([RotM*X'].');
   C0      = celestial.coo.cosined([RotM*X0'].');
   RA(J)   = C(1);
   Dec(J)  = C(2);
   RA0(J)  = C0(1);
   Dec0(J) = C0(2);
end
% The geocentric position angle of the northern semi minor axis of the
% apparent ellipse of the ring, measured from the north towerds
% the east.
P = atan((cos(Dec0).*sin(RA0-RA))./(sin(Dec0).*cos(Dec)-cos(Dec0).*sin(Dec).*cos(RA0-RA)));



