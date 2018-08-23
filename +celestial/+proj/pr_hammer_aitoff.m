function [X,Y]=pr_hammer_aitoff(Long,Lat,R)
%------------------------------------------------------------------------------
% pr_hammer_aitoff function                                           AstroMap
% Description: Project coordinates (longitude and latitude) using equal area
%              Hammer-Aitoff projection used in the FITS/WCS standard.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
%            For FITS WCS this parameter should be 180./pi.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                 September 2012     
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reference: Calabretta & Greisen 2002 A\&A 395, 1077
% Reliable: 1
%------------------------------------------------------------------------------
if (nargin==3),
   % no default
elseif (nargin==2),
   R = 1;
else
   error('Illigal number of argument');
end

% convert to [0 to pi] range:
Long = 2.*pi.*(Long./(2.*pi) - floor(Long./(2.*pi)));
% convert to [-pi to pi] range
I = find(Long>pi);
Long(I) = Long(I) - 2.*pi;

Gamma = R.*sqrt(2./(1+cos(Lat).*cos(Long.*0.5)));
X  = 2.*Gamma.*cos(Lat).*sin(0.5.*Long);
Y  = Gamma.*sin(Lat);

