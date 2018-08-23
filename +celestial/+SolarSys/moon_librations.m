function [LibMat]=moon_librations(JD,GeodCoo,MoonAppCoo)
% Moon's librations 
% Package: celestial.SolarSys
% Description: Calculate the Moon's librations and P.A., including
%                physical librations.
% Input  : - Vector of JD.
%          - Observer geodetic position, [East long, Lat, height],
%            in radians and meters.
%            NaN for geocentric librations.
%          - Optinal parameter : [Longitude, Latitude, HP] in radians.
%            where Longitude is a column vector of apparent geocentric
%            longitude of the Moon including the effect of the nutation.
%            where Latitude is a column vector of apparent geocentric
%            latitude of the Moon. Here, HP is a column vector of the
%            Moon's geocentric horiz. parallax. If not given it will
%            be calculated. If topocentric librations are needed then
%            topocentric coordinates should be given.
% Output : - [Libration in long, Libration in lat, P.A. of axis] in radians.
%            When the libration in longitude is positive, the the west limb
%            is exposed, when the libration in latitude is positive,
%            the north limb is exposed,
% Reference: Astro. Algo. J. Meeus 1991
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: BUGS
%------------------------------------------------------------------------------
RAD    = 180./pi;
FunTPI = inline(2.*pi.*(X./(2.*pi) - floor(X./(2.*pi)))','X');

if (nargin==1),
   GeodCoo    = NaN;
   [MoonRA,MoonDec,MoonHP] = celestial.SolarSys.mooncool(JD,GeodCoo,'b');

%   MoonAppCoo = 
if (nargin==2),
%   MoonAppCoo =  ... GeodCoo ...
elseif (nargin==3),
   % do nothing
else
   error('Illigal number of input arguments');
end

T = (JD - 2451545.0)./36525.0;
%
% The inclination of the mean lunar equator to the ecliptic (IAU value)
Inc = 1.54242./RAD;

% Moon apparent longitude
Long = MoonAppCoo(:,1);
% Moon apparent latitude
Lat  = MoonAppCoo(:,2);
% Moon H.P.
HP   = MoonAppCoo(:,3);

% calculate nutation in longitude (radians)
NutLon =


% argument of latitude of the Moon (radians)
F      = 93.2720993 + 483202.0175273.*T ...
        -0.0034029.*T.^2 - T.^3./3526000 + T.^4./863310000;
F      = F./RAD;
F      = FunTPI(F);

% mean longitude of ascending node of the lunar orbit (radians)
Omega  = 125.0445550 - 1934.1361849.*T ...
        +0.0020762.*T.^2 + T.^3./467410 - T.^4./60616000;
Omega  = Omega./RAD;
Omega  = FunTPI(Omega);

% Moon's mean anomaly
Mt     = 134.9634114 + 477198.8676313.*T ...
        +0.0089970.*T.^2 + T.^3./69699 - T.^4./14712000;
Mt     = Mt./RAD;
Mt     = FunTPI(Mt);

% Sun's mean anomaly
M      = 357.5291092 + 35999.0502909.*T ...
        -0.0001536.*T.^2 + T.^3./24490000;
M      = M./RAD;
M      = FunTPI(M);

% Earth eccentricity (change)
E      = 1 - 0.002516.*T - 0.0000074.*T.^2;

W = Long - NutLon - Omega;
TanA = (sin(W).*cos(Lat).*cos(Inc) - sin(Lat).*sin(Inc))./(cos(W).*cos(Lat));
A    = atan(TanA);
% optical libration in longitude
Ll   = A - F;

% optical libration in latitude
SinBl = -sin(W).*cos(Lat).*sin(I) - sin(Lat).*cos(I);
Bl    = asin(SinBl);


% calculate physical librations
K1    = (119.75 + 131.849.*T)./RAD;
K2    = ( 72.56 +  20.186.*T)./RAD;

Rho   = -0.02752.*cos(Mt) ...
        -0.02245.*sin(F) ...
        +0.00684.*cos(Mt - 2.*F) ...
        -0.00293.*cos(2.*F) ...
        -0.00085.*cos(2.*F - 2.*D) ...
        -0.00054.*cos(Mt - 2.*D) ...
        -0.00020.*sin(Mt + F) ...
        -0.00020.*cos(Mt + 2.*F) ...
        -0.00020.*cos(Mt - F) ...
        +0.00014.*cos(Mt + 2.*F - 2.*D);
Rho   = Rho./RAD;

Sigma = -0.02816.*sin(Mt) ...
        +0.02244.*cos(F) ...
        -0.00682.*sin(Mt - 2.*F) ...
        -0.00279.*sin(2.*F) ...
        -0.00083.*sin(2.*F - 2.*D) ...
        +0.00069.*sin(Mt - 2.*D) ...
        +0.00040.*cos(Mt + F) ...
        -0.00025.*sin(2.*Mt) ...
        -0.00023.*sin(Mt + 2.*F) ...
        +0.00020.*cos(Mt - F) ...
        +0.00019.*sin(Mt - F) ...
        +0.00013.*sin(Mt + 2.*F - 2.*D) ...
        -0.00010.*cos(Mt - 3.*F);
Sigma = Sigma./RAD;

Tau   = +0.02520.*sin(M).*E ...
        +0.00473.*sin(2.*Mt - 2.*F) ...
        -0.00467.*sin(Mt) ...
        +0.00396.*sin(K1) ...
        +0.00276.*sin(2.*Mt - 2.*D) ...
        +0.00196.*sin(Omega) ...
        -0.00183.*cos(Mt - F) ...
        +0.00115.*sin(Mt - 2.*D) ...
        -0.00096.*sin(Mt - D) ...
        +0.00046.*sin(2.*F - 2.*D) ...
        -0.00039.*sin(Mt - F) ...
        -0.00032.*sin(Mt - M - D) ...
        +0.00027.*sin(2.*Mt - M - 2.*D) ...
        +0.00023.*sin(K2) ...
        -0.00014.*sin(2.*D) ...
        +0.00014.*cos(2.*Mt - 2.*F) ...
        -0.00012.*sin(Mt - 2.*F) ...
        -0.00012.*sin(2.*Mt) ...
        +0.00011.*sin(2.*Mt - 2.*M -2.*D);
Tau   = Tau./RAD;

Phys_Ll = -Tau + (Rho.*cos(A) + Sigma.*sin(A)).*tan(Bl);
Phys_Bl = Sigma.*cos(A) - Rho.*sin(A);

TotLl = Ll + Phys_Ll;
TotBl = Bl + Phys_Bl; 

% Moon's position angle
Obl   = celestial.coo.obliquity(JD);
V     = Omega + NutLon + Sigma./sin(Inc);
X     = sin(Inc+Rho).*sin(V);
Y     = sin(Inc+Rho).*cos(V).*cos(Obl) - cos(Inc+Rho).*sin(Obl);
OmS   = atan2(X,Y);
Alpha =                 % moon's right ascen. 
SinP  = sqrt(X.^2 + Y.^2).*cos(Alpha - OmS)./cos(TotBl);
P     = asin(P);


LibMat = [TotLl, TotBl, P]; 




