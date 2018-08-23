function [X,Y]=pr_sinusoidal(Long,Lat,R)
%------------------------------------------------------------------------------
% pr_sinusoidal function                                              AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              Sinusoidal projection.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      July 1999      
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------
if (nargin==3)
   % no default
elseif (nargin==2)
   R = 1;
else
   error('Illigal number of argument');
end

% R is really R.*S (R-radius, S-scale factor)
X = Long.*R.*cos(Lat);
Y = R.*Lat;
