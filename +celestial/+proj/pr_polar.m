function [X,Y]=pr_polar(Long,Lat,R);
%------------------------------------------------------------------------------
% pr_polar function                                                   AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              polar projection (from north pole).
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    August 1999     
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin==3),
   % no default
elseif (nargin==2),
   % no default
   R = 1;
else
   error('Illigal number of argument');
end


CooLat = pi./2 - Lat;

X = R.*CooLat.*cos(Long);
Y = R.*CooLat.*sin(Long);



