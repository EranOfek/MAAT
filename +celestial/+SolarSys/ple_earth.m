function [L,B,R]=ple_earth(Date)
% Low-accuracy planetray ephemeris for Earth
% Package: celestial.SolarSys
% Description: Low accuracy planetray ephemeris for Earth. Calculate
%              Earth heliocentric longitude, latitude and radius
%              vector referred to the mean ecliptic and equinox of date.
%              Accuarcy: Better than 1' in long/lat, ~0.001 au in dist.
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
% Example: [L,B,R]=celestial.SolarSys.ple_earth([1 1 2000 0])
% Reliable: 2
%--------------------------------------------------------------------------


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

SumL0 = 175347046 ...
   + 3341656.*cos(4.6692568 + 6283.0758500.*Tau) ...
   + 34894.*cos(4.62610 + 12566.15170.*Tau) ...
   + 3497.*cos(2.7441 + 5753.3849.*Tau) ...
   + 3418.*cos(2.8289 + 3.5231.*Tau) ...
   + 3136.*cos(3.6277 + 77713.7715.*Tau) ...
   + 2676.*cos(4.4181 + 7860.4194.*Tau) ...
   + 2343.*cos(6.1352 + 3930.2097.*Tau) ...
   + 1324.*cos(0.7425 + 11506.7698.*Tau) ...
   + 1273.*cos(2.0371 + 529.6910.*Tau) ...
   + 1199.*cos(1.1096 + 1577.3435.*Tau) ...
   + 990.*cos(5.233 + 5884.927.*Tau) ...
   + 902.*cos(2.045 + 26.298.*Tau);

SumL1 = 628331966747 ...
   + 206059.*cos(2.678235 + 6283.075850.*Tau) ...
   + 4303.*cos(2.6351 + 12566.1517.*Tau) ...
   + 425.*cos(1.590 + 3.523.*Tau);

SumL2 = 52919 ...
   + 8720.*cos(1.0721 + 6283.0758.*Tau) ...
   + 309.*cos(0.867 + 12566.152.*Tau);

SumL3 = 289.*cos(5.844 + 6283.076.*Tau) ...
   + 35;

SumL4 = -114;

SumL5 = -1;


L = SumL0 + SumL1.*Tau + SumL2.*Tau.^2 ...
          + SumL3.*Tau.^3 + SumL4.*Tau.^4 + SumL5.*Tau.^5;
L = L.*1e-8;

L = FunTPI(L);

B = zeros(size(L));

%B = FunTPI(B);

SumR0 = 100013989 ...
   + 1670700.*cos(3.0984635 + 6283.0758500.*Tau) ...
   + 13956.*cos(3.05525 + 12566.15170.*Tau) ...
   + 3084.*cos(5.1985 + 77713.7715.*Tau) ...
   + 1628.*cos(1.1739 + 5753.3849.*Tau) ...
   + 1576.*cos(2.8469 + 7860.4194.*Tau);

SumR1 = 103019.*cos(1.107490 + 6283.075850.*Tau) ...
   + 1721.*cos(1.0644 + 12566.1517.*Tau) ...
   - 702;

SumR2 = 4359.*cos(5.7846 + 6283.0758.*Tau);

SumR3 = 145.*cos(4.273 + 6283.076.*Tau);

R = SumR0 + SumR1.*Tau + SumR2.*Tau.^2 ...
          + SumR3.*Tau.^3;
R = R.*1e-8;

