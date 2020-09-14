function [R,Vel,Nu]=elements2position(Date,OrbElem,Mu,GK)
% Convert orbital elements and date to position and velocity vectors.
% Package: celestial.Kepler
% Description: Convert orbital elements and date to position vector and
%              velocity vector at a give epoch.
% Input  : - Time of observations, one time per row.
%            If one column is given then taken as JD.
%            if four colum is given then taken as [D M Y frac_day].
%          - Orbital elements : [q e T i \Omega \omega], in au, days, radians.
%          - Gaussian gravitational constant.
%            Default is: 0.017202098950000
%          - GM [au^3 day^-2]. Default is: 2.959122084294439e-004
% Output : - Position matrix, [X;Y;Z] per column.
%            in au relative to ecliptic and equinox.
%          - Velocity matrix, [X;Y;Z] per column.
%            in au/day relative to ecliptic and equinox.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.Kepler.elements2position(2451545,[1 0.1 2451545 0 0 0])
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==2)
   % Gaussian grav. constant
   GK = 0.017202098950000;
   Mu = 2.959122084294439e-004;   
elseif (nargin==3)
   %AU = 1.49597870660e11; % meters
   Mu = 2.959122084294439e-004;
   %Mu = 1.32712440e20./(AU.^3)./(86400.^2);
elseif (nargin==4)
   % do nothing
else
   error('Illigal number of input arguments');
end

SizeDate = size(Date);
Len      = SizeDate(1);
ColN     = SizeDate(2);
if (ColN==1)
   % already in JD
   JD = Date;
elseif (ColN==4)
   JD = julday(Date);
else
   error('Illigal number of columns in Date matrix');
end

Q     = OrbElem(1);
E     = OrbElem(2);
T     = OrbElem(3);
Inc   = OrbElem(4);
Omega = OrbElem(5);
W     = OrbElem(6);

if (E<1 && E>=0)
   % elliptic motion
   %A      = Q.*(1 - E);
   %MM  = sqrt(Mu./(A.^3));   % mean motion [rad/day]
   %M   = MM.*(JD - T);      % Mean anomaly
   [Nu,Rad,~,Vel] = celestial.Kepler.kepler_elliptic(JD-T,Q,E,GK);
   
elseif (E==1)
   % parabolic motion
   [Nu,Rad,~,Vel] = celestial.Kepler.kepler_parabolic(JD-T,Q,GK);
   
elseif (E>1)
   % hyperbolic equation
   [Nu,Rad,~,Vel] = celestial.Kepler.kepler_hyperbolic(JD-T,Q,E,GK);
  
else
   error('Illigal eccentricity, e<0');
end

[X,Y,Z] = celestial.Kepler.trueanom2pos(Rad,Nu,Omega,W,Inc);

R = [X,Y,Z];
