function varargout=trueanom2pos(R,Ni,Omega,OmP,Inc)
% True anomaly, radius vector and orbital elements to cartezian position
% Package: celestial.Kepler
% Description: Given an object true anomaly, radius vector, time and
%              orbital elements and time, calculate its orbital position
%              in respect to the orbital elements reference frame.
% Input  : - Vector of radius vector, [length unit].
%          - Vector of True anomaly, [radians].
%          - Argument of Ascending node, [radians].
%          - Longitude of periastron, [radians]. 
%          - Inclination, [radians].
% Output : * If one output argument is requested then the program returns
%            a single matrix containing the [X,Y,Z].
%            If three output arguments are requested them the program
%            returns the X, Y, Z position vectors in ecliptic coordinate
%            system.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: trueanom2vel.m
% Example: [X,Y,Z] = celestial.Kepler.trueanom2pos(1,[0;1],1,1,1);
%          [XYZ] = celestial.Kepler.trueanom2pos(1,[0;1],1,1,1);
%          [XYZ] = celestial.Kepler.trueanom2pos([1;2],[0;1],1,1,1);
%          [XYZ] = celestial.Kepler.trueanom2pos([1;2],[0;1],1,1,[1;0.1]);
% Reliable: 1
%------------------------------------------------------------------------------

X = R.*(cos(OmP+Ni).*cos(Omega) - sin(OmP+Ni).*sin(Omega).*cos(Inc));
Y = R.*(cos(OmP+Ni).*sin(Omega) + sin(OmP+Ni).*cos(Omega).*cos(Inc));
Z = R.*sin(OmP+Ni).*sin(Inc);

if (nargout==1)
   varargout{1} = [X,Y,Z];
else
   varargout{1} = X;
   varargout{2} = Y;
   varargout{3} = Z;
end
