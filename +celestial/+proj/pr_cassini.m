function [X,Y]=pr_cassini(Long,Lat,R,CenLong)
% Cassini projection.
% Package: celestial.proj
% Description: Project coordinates (longitude and latitude) using the
%              Cassini projection.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
%          - Central longitude, defauly=t is 0.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jul 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

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


B = cos(Lat).*sin(Long - CenLong);
X = R.*asin(B);
Y = R.*atan(tan(Lat)./cos(Long-CenLong));
