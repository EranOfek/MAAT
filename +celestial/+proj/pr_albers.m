function [X,Y]=pr_albers(Long,Lat,R,CenCoo,ParLat)
% Albers Equal-Area projection.
% Package: celestial.proj
% Description: Project coordinates (longitude and latitude) using the
%              Albers Equal-Area projection. The coordinates are projected on
%              ellipse with axis ratio of 2:1.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
%          - The origin point of the Cartesian coordinates
%            [Long_0, Lat_0], default is [0 0].
%          - The standart parallels [Lat_1, Lat_2],
%            default is [-pi/2 pi/2].
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jul 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
if (nargin==2),
   R = 1;
   CenCoo = [0 0];
   ParLat = [-pi/2 pi/2];
elseif (nargin==3),
   CenCoo = [0 0];
   ParLat = [-pi/2 pi/2];
elseif (nargin==4),
   ParLat = [-pi/2 pi/2];
elseif (nargin==5),
   % no default   
else
   error('Illigal number of argument');
end

Lat_1  = ParLat(1);
Lat_2  = ParLat(2);
Long_0 = CenCoo(1);
Lat_0  = CenCoo(2);

n     = 0.5.*(sin(Lat_1) + sin(Lat_2));
C     = cos(Lat_1).^2 + 2.*n.*sin(Lat_0);
Rho_0 = sqrt(C-2.*n.*sin(Lat_0));
Theta = n.*(Long-Long_0);
Rho   = sqrt(C-2.*n.*sin(Lat))./n;

X = Rho.*sin(Theta);
Y = Rho_0 - Rho.*cos(Theta);

X = R.*X;
Y = R.*Y;
