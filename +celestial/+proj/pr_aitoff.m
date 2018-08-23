function [X,Y]=pr_aitoff(Long,Lat,R)
% Project coordinates using equal area Aitoff projection
% Package: celestial.proj
% Description: Project coordinates (longitude and latitude) using equal
%              area Aitoff projection.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 1999     
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [X,Y]=celestial.proj.pr_aitoff(1,1,1)
% Reliable: 1
%--------------------------------------------------------------------------
if (nargin==3),
   % no default
elseif (nargin==2),
   R = 1;
else
   error('Illigal number of argument');
end

% R is really R.*S (R-radius, S-scale factor)
%X = -2.*R.*cos(Lat).*sin(0.5.*Long)./sqrt(1+cos(Lat).*cos(0.5.*Long));
%Y = R.*sin(Lat)./sqrt(1+cos(Lat).*cos(0.5.*Long));

Alpha = acos(cos(Lat).*cos(0.5.*Long));
X  = 2.*R.*cos(Lat).*sin(0.5.*Long)./sinc(Alpha./pi);
Y  = sin(Lat)./sinc(Alpha./pi);
