function [X,Y]=pr_gnomonic(Long,Lat,R,CenterVec)
%--------------------------------------------------------------------------
% pr_gnomonic function                                            AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              Gnomonic non conformal projection
%              This is a nonconformal projection from a sphere center in
%              which orthodromes are stright lines.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
%          - Central coordinate vector [Long_center,Lat_center]
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jul 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: pr_ignomonic.m
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==4)
   % no default
elseif (nargin==3)
   CenterVec = [0 0];
elseif (nargin==2)
   CenterVec = [0 0];
   R = 1;
else
   error('Illigal number of argument');
end
%CenterVec = [-0.2514 0.4727];
% extract central coordinates.

% Change CenterVec by data middle
%--------------------------------
%CenterVec(1) = (max(Long)+min(Long)).*0.5;
%CenterVec(2) = (max(Lat)+min(Lat)).*0.5;
%--------------------------------

Long1 = CenterVec(1);
Lat1  = CenterVec(2);
% R is really R.*S (R-radius, S-scale factor)
CosC = sin(Lat1).*sin(Lat) + cos(Lat1).*cos(Lat).*cos(Long-Long1);
X = R.*cos(Lat).*sin(Long-Long1)./CosC;
Y = R.*(cos(Lat1).*sin(Lat) - sin(Lat1).*cos(Lat).*cos(Long-Long1))./CosC;

