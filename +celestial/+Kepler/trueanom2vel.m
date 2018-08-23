function [X,Y,Z,Xd,Yd,Zd]=trueanom2vel(R,Ni,Rd,Nid,Omega,OmP,Inc)
% True anomaly, radius vector and orbital elements to position and velocity
% Package: celestial.Kepler
% Description: Given an object true anomaly, radius vector, their
%              derivatives and orbital elements and time, calculate its
%              orbital position and velocity in respect to the orbital
%              elements reference frame.
% Input  : - Vector of radius vector, [length unit].
%          - Vector of True anomaly, [radians].
%          - Vector of radius vector derivative [radians/time].
%          - Vector of true anomaly derivative [radians/time].
%          - Argument of Ascending node, [radians].
%          - Longitude of periastron, [radians]. 
%          - Inclination, [radians].
% Output : - X
%          - Y
%          - Z
%          - X dot
%          - Y dot
%          - Z dot
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: trueanom2pos.m
% Example: [X,Y,Z,Xd,Yd,Zd] = celestial.Kepler.trueanom2vel(1,[0;1],0.1,0.1,1,1,1);
% Reliable: 2
%------------------------------------------------------------------------------
'require testing!!' 
X = R.*(cos(OmP+Ni).*cos(Omega) - sin(OmP+Ni).*sin(Omega).*cos(Inc));
Y = R.*(cos(OmP+Ni).*sin(Omega) + sin(OmP+Ni).*cos(Omega).*cos(Inc));
Z = R.*sin(OmP+Ni).*sin(Inc);

Xd = Rd.*X./R - R.*Nid.*(sin(OmP+Ni).*cos(Omega) + cos(OmP+Ni).*sin(Omega).*cos(Inc));
Yd = Rd.*Y./R - R.*Nid.*(sin(OmP+Ni).*sin(Omega) - cos(OmP+Ni).*cos(Omega).*cos(Inc));
Zd = Rd.*Z./R + R.*Nid.*cos(OmP+Ni).*sin(Inc);


