function [Lon,Lat]=celestial_circ(CenLon,CenLat,Rad,Np)
% Grid of coordinates on a small spherical circle
% Package: celestial.coo
% Description: Calculate grid of longitude and latitude of a small circle
%              on the celestial sphere.
% Input  : - Longitude of small circle center [radians].
%          - Latitude of small circle center [radians].
%          - Radius of small circle [radians].
%          - Number of points to calculate, default is 100.
% Output : - Longitude of points on the small circle.
%          - Latitude of points on the small circle.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Dec 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Lon,Lat]=celestial.coo.celestial_circ(1,1,1,500);
% Reliable: 2
%--------------------------------------------------------------------------
if (nargin==3)
   Np = 100;
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

Az         = 2.*pi.*(1:1:Np).'./Np;
Ones       = ones(Np,1);
[Lat, Lon] = reckon(CenLat.*Ones,CenLon.*Ones,Rad.*Ones,Az,'radians');
