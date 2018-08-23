function [L,B,R]=ple_venus(Date)
% Low accuracy ephemeris for Venus
% Package: celestial.SolarSys
% Description: Low accuracy ephemeris for Venus. Calculate 
%              Venus heliocentric longitude, latitude and radius 
%              vector referred to mean ecliptic and equinox of date.
%              Accuarcy:  Better than 1' in long/lat, ~0.001 au in dist.
% Input  : - matrix of dates, [D M Y frac_day] per line,
%            or JD per line. In TT time scale.
% Output : - Longitude in radians.
%          - Latitude in radians.
%          - Radius vector in au.
% Reference: VSOP87
% See also: ple_planet.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [L,B,R]=celestial.SolarSys.ple_venus(2451545)
% Reliable: 2
%------------------------------------------------------------------------------


RAD = 180./pi;

FunTPI = @(X) (X./(2.*pi) - floor(X./(2.*pi))).*2.*pi;

SizeDate = size(Date);
N        = SizeDate(1);
ColN     = SizeDate(2);

if (ColN==4),
   JD = celestial.time.julday(Date);
elseif (ColN==1),
   JD = Date;
else
   error('Illigal number of columns in date matrix');
end

Tau = (JD(:) - 2451545.0)./365250.0;

SumL0 = 317614667 ...
   + 1353968.*cos(5.5931332 + 10213.2855462.*Tau) ...
   + 89892.*cos(5.30650 + 20426.57109.*Tau) ...
   + 5477.*cos(4.4163 + 7860.4194.*Tau) ...
   + 3456.*cos(2.6996 + 11790.6291.*Tau) ...
   + 2372.*cos(2.9938 + 3930.2097.*Tau) ...
   + 1664.*cos(4.2502 + 1577.3435.*Tau) ...
   + 1438.*cos(4.1575 + 9683.5946.*Tau) ...
   + 1317.*cos(5.1867 + 26.2983.*Tau) ...
   + 1201.*cos(6.1536 + 30639.8566.*Tau) ...
   + 769.*cos(0.816 + 9437.763.*Tau);

SumL1 = 1021352943053 ...
   + 95708.*cos(2.46424 + 10213.28555.*Tau) ...
   + 14445.*cos(0.51625 + 20426.57109.*Tau);

SumL2 = 54127 ...
   + 3891.*cos(0.3451 + 10213.2855.*Tau) ...
   + 1338.*cos(2.0201 + 20426.5711.*Tau);

SumL3 = 136.*cos(4.804 + 10213.286.*Tau);

SumL4 = -114;

SumL5 = -1;


L = SumL0 + SumL1.*Tau + SumL2.*Tau.^2 ...
          + SumL3.*Tau.^3 + SumL4.*Tau.^4 + SumL5.*Tau.^5;
L = L.*1e-8;

L = FunTPI(L);

SumB0 = 5923638.*cos(0.2670278 + 10213.2855426.*Tau) ...
   + 40108.*cos(1.14737 + 20426.57109.*Tau) ...
   - 32815 ...
   + 1011.*cos(1.0895 + 30639.8566.*Tau);

SumB1 = 513348.*cos(1.803643 + 10213.285546.*Tau) ...
   + 4380.*cos(3.3862 + 20426.5711.*Tau) ...
   + 199;

SumB2 = 22378.*cos(3.38509 + 10213.28555.*Tau) ...
   + 282;

SumB3 = 647.*cos(4.992 + 10213.286.*Tau) ...
   - 20;

SumB4 = 44.*cos(0.32 + 10213.29.*Tau);

B = SumB0 + SumB1.*Tau + SumB2.*Tau.^2 ...
          + SumB3.*Tau.^3 + SumB4.*Tau.^4;
B = B.*1e-8;

%B = FunTPI(B);

SumR0 = 72334821 ...
   + 489824.*cos(4.021518 + 10213.285546.*Tau) ...
   + 1658.*cos(4.9021 + 20426.5711.*Tau);

SumR1 = 34551.*cos(0.89199 + 10213.28555.*Tau);

SumR2 = 1407.*cos(5.0637 + 10213.2855.*Tau);

SumR3 = 50.*cos(3.22 + 10213.29.*Tau);

R = SumR0 + SumR1.*Tau + SumR2.*Tau.^2 ...
          + SumR3.*Tau.^3;
R = R.*1e-8;

