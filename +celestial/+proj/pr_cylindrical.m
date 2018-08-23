function [X,Y]=pr_cylindrical(Long,Lat,R,StandCoo);
%------------------------------------------------------------------------------
% pr_cylindrical function                                             AstroMap
% Description: Project coordinates (longitude and latitude) using a general
%              cylindrical projection.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
%          - Standard coordinates. [Stand_Long, Stand_Lat]. (radians)
%            Speicial cases of cylindrical equal-area projections
%            Stand_Lat    Projection name
%            0 deg.       Lambert Cylindrical Equal-Area (default)
%            30 deg.      Behrmann Cylindrical Equal-Area
%            37.383 deg.  Tristan Edwards Projection
%            44.138 deg.  Peters Projection
%            45 deg.      Gall Orthographic projection
%            50 deg.      Balthasart projection
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      July 1999
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------
if (nargin==2),
   R = 1;
   StandCoo = [0 0];
elseif (nargin==3),
   StandCoo = [0 0];
elseif (nargin==4),
   % no default   
else
   error('Illigal number of argument');
end

StandLong = StandCoo(1);
StandLat  = StandCoo(2);

X = R.*(Long - StandLong).*cos(StandLat);
Y = R.*sin(Lat).*sec(StandLat);
