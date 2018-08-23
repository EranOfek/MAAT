function [RA,Dec,RadVec]=suncoo1(JD)
% Low-accuracy coordinates of the Sun (1950-2050 range)
% Package: celestial.SolarSys
% Description: Calculate the Sun equatorial coordinates using
%              low accuracy formaulae for the range 1950 to 2050.
%              Accuracy : 0.01 deg. in long, 0.1m in Equation of Time
% Input  : - Vector of JDs
% Output : - vector of RA, in radians.
%          - vector of Dec. in radians.
%          - Vector of radius vectors, in AU.
% Tested : Matlab 5.2
%     By : Eran O. Ofek                    August 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: suncoo.m, mooncoo.m
% Example: [RA,Dec,RadVec]=celestial.SolarSys.suncoo1(julday([1 1 2000]));
% Reliable: 1
%------------------------------------------------------------------------------
RAD = 180./pi;

n  = JD - 2451545.0;
if (max(abs(n))>(50.*365)),
   error('This formulae give good results only between 1950-2050');
end

L   = (280.466 + 0.9856474.*n)./RAD;
g   = (357.528 + 0.9856003.*n)./RAD;

Lam = L + (1.915.*sin(g) + 0.020.*sin(2.*g))./RAD;
Bet = 0;
Obl = (23.440 - 0.0000004.*n)./RAD;

RadVec = 1.00014 - 0.01671.*cos(g) - 0.00014.*cos(2.*g);
RA     = atan2(cos(Obl).*sin(Lam),cos(Lam));
Dec    = asin(sin(Obl).*sin(Lam));
