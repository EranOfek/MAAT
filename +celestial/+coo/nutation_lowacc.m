function [DLon,DObl]=nutation_lowacc(JD)
%--------------------------------------------------------------------------
% nutation_lowacc function                                           ephem
% Description: Low accuracy (~1") calculation of the nutation.
% Input  : - Matrix of Julian days.
% Output : - Matrix of nutation in longitude [rad].
%          - Matrix of nutation in obliquity [rad].
% See also: nutation.m, nutation1984.m, nutation2rotmat.m
% Reference: Seidelmann (1992)
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DLon,DObl]=nutation_lowacc(2451666);
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;
D = JD - 2451545.0;

DLon = (-0.0048.*sind(125.0-0.05295.*D) - 0.0004.*sind(200.9+1.97129.*D))./RAD;
DObl = (+0.0026.*cosd(125.0-0.05295.*D) + 0.0002.*cosd(200.9+1.97129.*D))./RAD;

