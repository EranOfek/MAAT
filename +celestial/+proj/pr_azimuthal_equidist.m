function [X,Y]=pr_azimuthal_equidist(Long,CoLat,R)
% Azimuthal equidistant projection.
% Package: celestial.proj
% Description: Project longitude and latitude using Azimuthal equidistant
%              projection (constant radial scale).
% Input  : - Longitude [rad].
%          - CoLatitude [rad].
%          - Sphere radius, default is 1.
% Output : - X
%          - Y
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Aug 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%------------------------------------------------------------------------------

RAD = 180./pi;

if (nargin==2),
   R = 1;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end


Radius = R.*CoLat./(pi./2);

X = Radius.*cos(Long);
Y = Radius.*sin(Long);
