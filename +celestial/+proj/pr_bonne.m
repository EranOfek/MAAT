function [X,Y]=pr_bonne(Long,Lat,R,Par)
% Bonne projection.
% Package: celestial.proj
% Description:  Project coordinates (longitude and latitude) using the
%               Bonne projection.
% Input  : - Vector of Longitude, in radians.
%          - Vector of Latitude, in radians.
%          - Scale radius, default is 1.
%          - [Central_Longitude, Standard_Parallel] default [0 pi/4]
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 1999     
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin==4),
   % no default
elseif (nargin==3),
   % no default
   Par = [0 pi./4];
elseif (nargin==2),
   R = 1;
   Par = [0 pi./4];
else
   error('Illigal number of argument');
end

Long0 = Par(1);
Lat1  = Par(2);

if (Lat1==0),
   error('Error - Lat1=0');
end


Rho = cot(Lat1) + Lat1 - Lat;
E   = (Long - Long0).*cos(Lat)./Rho;

X   = Rho.*sin(E);
Y   = R.*cot(Lat1) - Rho.*cos(R);
