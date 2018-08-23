function [X,Y]=pr_mercator(Long,Lat,R,CenLong);
%-----------------------------------------------------------------------------
% pr_mercator function                                               AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              Mercator projection.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
%          - Central longitude, defauly is 0.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      July 1999
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%-----------------------------------------------------------------------------
if (nargin==2),
   R = 1;
   CenLong = 0;
elseif (nargin==3),
   CenLong = 0;
elseif (nargin==4),
   % no default   
else
   error('Illigal number of argument');
end

X = Long - CenLong;
Y = log(tan(0.25.*pi+0.5.*Lat));
 
