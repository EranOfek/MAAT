function [RA,Dec,R,SL,EquationTime]=suncoo(JD,EquinoxType)
% Low-accuracy position of the Sun (0.01 deg in long).
% Package: celestial.SolarSys
% Description: Calculate the Sun equatorial coordinates using low
%              accuracy formale. Accuracy : 0.01 deg. in long.
% Input  : - Vector of JDs.
%          - Equinox for output coordinates:
%            'g' - True Geometric referred to the mean equinox of date.
%            'a' - Apparent, referred to the true equinox of date. (default).
%            'j' - J2000, referred to the J2000.0 equinox.
% Output : - vector of RA, in radians.
%          - vector of Dec. in radians.
%          - Vector of radius vectors, in AU.
%          - Solar longitude in the same ref. frame as RA/Dec. (radians).
%          - Equation of Time [Minuts of time]
% See Also : suncoo1; mooncoo
% Tested : matlab 5.3
%     By : Eran O. Ofek                    Sep 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [RA,Dec,R,SL,EquationTime]=celestial.SolarSys.suncoo(2451545+[0:1:10]','j');
% Reliable: 1
%------------------------------------------------------------------------------
if (nargin==1),
   EquinoxType = 'a';
elseif (nargin==2),
   % no default
else
   error('Illigal number of input arguments');
end
RAD = 180./pi;

T   = (JD - 2451545.0)./36525;
L0  = (280.46645 + 36000.76983.*T + 0.0003032.*T.*T)./RAD;
M   = (357.52910 + 35999.05030.*T - 0.0001559.*T.*T - 0.00000048.*T.*T.*T)./RAD;
e   = 0.016708617 - 0.000042037.*T - 0.0000001236.*T.*T;

C   = (1.914600 - 0.004817.*T - 0.000014.*T.*T).*sin(M) + (0.019993 - 0.000101.*T).*sin(2.*M) + 0.000290.*sin(3.*M);
C   = C./RAD;
% Sun longitude
SL  = L0 + C;

% the sun true Anomaly:
Ni  = M + C;

% the sun radius vector
R   = 1.000001018.*(1-e.^2)./(1+e.*cos(Ni));

if (EquinoxType=='a'),
   Om = (125.04 - 1934.136.*T)./RAD;
   SL = SL - (0.00569 - 0.00478.*sin(Om))./RAD;
elseif (EquinoxType=='j'),
   SL = SL - (0.01397.*T.*100)./RAD;
elseif (EquinoxType=='g'),
   % Allready geometric longitude
else
   error('Illigal equinox type');
end


Obl = celestial.coo.obliquity(JD);

SL     = (SL./(2.*pi) - floor(SL./(2.*pi))).*2.*pi;
RA     = atan2(cos(Obl).*sin(SL),cos(SL));
Dec    = asin(sin(Obl).*sin(SL));


EquationTime = L0 - 0.0057183./RAD - RA;
EquationTime = 1440.*(EquationTime./(2.*pi) - floor(EquationTime./(2.*pi)));
if (EquationTime>720),
   EquationTime = EquationTime - 1440;
end

