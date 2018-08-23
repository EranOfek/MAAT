function [X,Y]=pr_xy(Long,Lat,R)
% X-Y projection (no transformation).
% Package: celestial.proj
% Description: Project coordinates (longitude and latitude) to X-Y
%              projection (no transformation).
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.2
%     By : Eran O. Ofek                    Aug 1999     
%    URL : http://weizmann.ac.il/home/eofek/matlab/
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


X = R.*Long./(2.*pi);
Y = R.*Lat./(pi);



