function Area=spherical_tri_area(Coo1,Coo2,Coo3);
%------------------------------------------------------------------------------
% spherical_tri_area function                                         AstroMap
% Description: Given three coordinates on a sphere, calculate the area
%              of a spherical triangle defined by these three points.
% Input  : - [Long, Lat] in radians of first point on a sphere.
%          - [Long, Lat] in radians of second point on a sphere.
%          - [Long, Lat] in radians of third point on a sphere.
% Output : - Area of spherical triangles.
% Tested : Matlab 7.6
%     By : Eran O. Ofek                  December 2008
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

Cos.Coo1 = cosined(Coo1);
Cos.Coo2 = cosined(Coo2);
Cos.Coo3 = cosined(Coo3);

% angles between points
Ang12 = acos(dot(Cos.Coo1,Cos.Coo2,2));  
Ang23 = acos(dot(Cos.Coo2,Cos.Coo3,2));
Ang31 = acos(dot(Cos.Coo3,Cos.Coo1,2));

% opposite angles
Opp12 = acos((cos(Ang12) - cos(Ang23).*cos(Ang31))./(sin(Ang23).*sin(Ang31)));
Opp23 = acos((cos(Ang23) - cos(Ang31).*cos(Ang12))./(sin(Ang31).*sin(Ang12)));
Opp31 = acos((cos(Ang31) - cos(Ang12).*cos(Ang23))./(sin(Ang12).*sin(Ang23)));

% area
Area  = Opp12 + Opp23 + Opp31 - pi;
