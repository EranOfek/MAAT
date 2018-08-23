function [X,Y]=pr_hammer(Long,Lat,R);
%------------------------------------------------------------------------------
% pr_hammer function                                                  AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              Hammer projection. The coordinates are projected on an ellipse
%              with axis ratio of 2:1.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    August 1999     
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 1
%------------------------------------------------------------------------------
if (nargin==3),
   % no default
elseif (nargin==2),
   R = 1;
else
   error('Illigal number of argument');
end

X = 2.*R.*sqrt(2).*cos(Lat).*sin(Long./2)./sqrt(1+cos(Lat).*cos(Long./2));
Y = R.*sqrt(2).*sin(Lat)./sqrt(1+cos(Lat).*cos(Long./2));
