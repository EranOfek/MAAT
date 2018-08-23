function [X,Y]=pr_parabolic(Long,Lat,R);
%------------------------------------------------------------------------------
% pr_parabolic function                                               AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              Parabolic projection.
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
if (nargin==3),
   % no default
elseif (nargin==2),
   R = 1;
else
   error('Illigal number of argument');
end

% R is really R.*S (R-radius, S-scale factor)
X = 1.53499.*2.*Long.*R.*(2.*cos(2.*Lat./3)-1)./pi;
Y = 3.06998.*R.*sin(Lat./3);
