function [X,Y]=pr_mollweide(Long,Lat,R);
%------------------------------------------------------------------------------
% pr_mollweide function                                               AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              equal area Mollweide projection.
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

% R is really R.*S (R-radius, S-scale factor)
X = 2.*Long.*sqrt(2).*R.*cos(Lat)./pi;
Y = sqrt(2).*R.*sin(Lat);
