function [X,Y]=pr_conic(Long,Lat,R,Par)
% Conic projection.
% Package: celestial.proj
% Description: Project coordinates (longitude and latitude) using the
%              Conic projection.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
%          - Height of apex above sphare center, in sphare radius units.
%            default is 2.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 1999     
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==4),
   % no default
elseif (nargin==3),
   % no default
   Par = 2;
elseif (nargin==2),
   R = 1;
   Par = 2;
else
   error('Illigal number of argument');
end


H  = Par;
%SH = sqrt(H.*H - 1);
%X = csc(asec(H) + Lat).*cos(Lat).*sin(Long./SH);
%Y = csc(asec(H) + Lat).*cos(Lat).*cos(Long./SH);

CooLat = pi./2 - Lat;
Theta  = asin(R./H);

L   = H.*sin(Theta)./sin(pi - Theta - CooLat);
Rad = L.*sin(pi-CooLat);
%Rad2 = L.*L + H.*H - 2.*L.*H.*cos(CooLat);
%Rad  = sqrt(Rad2);

X   = Rad.*cos(Long);
Y   = Rad.*sin(Long);
