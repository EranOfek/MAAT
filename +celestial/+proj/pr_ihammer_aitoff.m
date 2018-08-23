function [Long,Lat]=pr_ihammer_aitoff(X,Y,R)
%------------------------------------------------------------------------------
% pr_ihammer_aitoff function                                          AstroMap
% Description: Project coordinates (longitude and latitude) using the inverse
%              of the equal area Hammer-Aitoff projection used in the
%              FITS/WCS standard.
% Input  : - Vector of X position.
%          - Vector of X position.
%          - Scale radius, default is 1.
%            For FITS WCS this parameter should be 180./pi.
% Output : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                 September 2012     
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


Z = sqrt(1 - (X./(4.*R)).^2 - (Y./(2.*R)).^2);

Lat  = asin(Y.*Z./R);
Long = 2.*atan2(Z.*X./(2.*R), 2.*Z.^2-1);

