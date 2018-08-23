function [X,Y]=pr_planis(Long,Lat,Phi,R)
% planisphere projection.
% Package: celestial.proj
% Description: Project longitude and latitude using a
%              'planisphere projection'.
% Input  : - Longitude [rad].
%          - Latitude [rad].
%          - Latitude [rad] defining projection point, default is 32./RAD.
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
   Phi = 32./RAD;
   R = 1;
elseif (nargin==3),
   R = 1;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end



Z = pi./2 - Lat;
L = R + R./sin(pi./2 - Phi);
T = sqrt(R.^2 + (L-R).^2 - 2.*R.*(L-R).*cos(pi - Z));
Alpha = asin(R.*sin(pi - Z)./T);
Radius = L.*tan(Alpha);

X = Radius.*cos(Long);
Y = Radius.*sin(Long);
