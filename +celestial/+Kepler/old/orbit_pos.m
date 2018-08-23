function [X,Y,Z]=orbit_pos(t,Elem,Mass2,n)
%--------------------------------------------------------------------------
% orbit_pos function                                                 ephem
% Description: Given an object true anomaly, radius vector, time and
%              orbital elements and time, calculate its orbital position
%              in respect to the orbital elements reference frame.
% Input  : - Vector of time [days].
%          - Orbital elements : [T, q, e, Om, w, i],
%            The angles are in radians and time in days.
%          - The mass of the secondary in units of the
%            primary mass's.
%            default is the earth+moon mass 1./328900.5 
%          - Optional mean motion - if given then mass is not used.
% Output : - X, rectangular position.
%          - Y
%          - Z
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    August 1999
%    URL : http://wise-obs.tau.ac.il~/eran/matlab.html
%--------------------------------------------------------------------------

if (nargin==2),
   Mass2 = 1./328900.5;
   n     = [];
elseif (nargin==3),
   n     = [];
elseif (nargin==4),
   % default
else
   error('Illigal number of arguments');
end

T  = Elem(:,1);
q  = Elem(:,2);
e  = Elem(:,3);
Om = Elem(:,4);
w  = Elem(:,5);
i  = Elem(:,6);

X  = zeros(length(t),1);
Y  = zeros(length(t),1);
Z  = zeros(length(t),1);


for I=1:1:length(t),
   [r,ni] = kepler(t,Elem,Mass2,n);
   % rectangular coordinates
   % in case of Solar system, heliocentric referred to the ecliptic of date.
   X      = r.*(cos(ni+w).*cos(Om) - sin(ni+w).*cos(i).*sin(Om));    
   Y      = r.*(cos(ni+w).*sin(Om) + sin(ni+w).*cos(i).*cos(Om));
   Z      = r.*sin(ni+w).*sin(i);
end
