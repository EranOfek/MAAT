function [Elem,Elem2]=position2elements(Date,R,V,Mu,GK)
% Convert position/velocity vectors to orbital elements.
% Package: celestial.Kepler
% Description: Convert date, position vector and velocity vector at a give
%              ephoch to orbital elements at the same epoch.
% Input  : - Time of observations, one time per raw.
%            If one column is given then taken as JD.
%            if four colum is given then taken as [D M Y frac_day].
%          - Position matrix, [X;Y;Z] per column.
%            in au relative to ecliptic and equinox.
%          - Velocity matrix, [X;Y;Z] per column.
%            in au/day relative to ecliptic and equinox.
%          - Gaussian gravitational constant.
%            Default is: 0.017202098950000
%          - GM [au^3 day^-2]. Default is: 2.959122084294439e-004
% Output : - [q, e, T, i, \Omega, \omega], distance in au, time in days,
%            angles in radians.
%          - [a, \nu, EccentricAnomaly, MeanAnomaly, MeanLongitude, n], distance
%            in au, time in days, angles in radians.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: NOT READY
%------------------------------------------------------------------------------

error('Not ready for use');

if (nargin==3),
   % Gaussian grav. constant
   GK = 0.017202098950000;
   Mu = 2.959122084294439e-004;   
elseif (nargin==4),
   AU = 1.49597870660e11; % meters
   Mu = 2.959122084294439e-004;
   %Mu = 1.32712440e20./(AU.^3)./(86400.^2);
elseif (nargin==5),
   % do nothing
else
   error('Illigal number of input arguments');
end

SizeDate = size(Date);
Len      = SizeDate(1);
ColN     = SizeDate(2);
if (ColN==1),
   % already in JD
   JD = Date;
elseif (ColN==4),
   JD = julday(Date);
else
   error('Illigal number of columns in Date matrix');
end

Dim = 1;
AbsVal = inline('sqrt(sum(X.*X,Dim))','X','Dim');
AbsR = AbsVal(R,Dim);
AbsV = AbsVal(V,Dim);
RV   = sum(R.*V,Dim);

% eccentricity
E    = (AbsV.^2./Mu - 1./AbsR).*R - (RV./Mu).*V;
AbsE = AbsVal(E,Dim);

% Angular momentum vector
H    = cross(R,V,Dim);
AbsH = AbsVal(H,Dim);

% Ascending node vector
if (Dim==1),
   K = [0;0;1]*ones(1,Len);
   I = [1;0;0]*ones(1,Len);   
else
   K = ones(Len,1)*[0 0 1];
   I = ones(Len,1)*[1 0 0];   
end
AbsK = AbsVal(K,Dim);
AbsI = AbsVal(I,Dim);
N    = cross(K,H);
AbsN = AbsVal(N,Dim);

% semi major axis
A   = 1./(2./AbsR - AbsV.^2./Mu);
% semi latus
LR  = AbsH.^2./Mu;

% periastron distance
Q   = LR./(1 + AbsE);

Inc   = acos(sum(K.*H,Dim)./(AbsK.*AbsH));
Omega = acos(sum(I.*N,Dim)./(AbsI.*AbsN));
INy   = find(N(2,:)<0);
Omega(INy) = 2.*pi - Omega(INy);
W     = acos(sum(N.*E,Dim)./(AbsN.*AbsE));
IeZ   = find(E(3,:)<0);
W(IeZ)= 2.*pi - W(IeZ);

B     = A.*sqrt(1 - AbsE.^2);

MeanX = (LR - AbsR)./AbsE;
MeanY = (RV./AbsE).*sqrt(LR./Mu);
% True anomaly
Nu    = atan2(MeanY./AbsR, MeanX./AbsR);
% Eccentric anomaly
CosE = MeanX./A + AbsE; 
SinE = MeanY./B;
EccA = atan2(SinE, CosE);
% Mean anomaly
if (AbsE>1),
   % hyperbolic (e>1)
   B    = -A.*sqrt(AbsE.^2 - 1);
   SinH = MeanY./B;
   M    = AbsE.*SinH - asin(SinH);
   % Mean motion
   MM   = GK.*sqrt(Mu./(-A.^3));
elseif (abs(AbsE-1)<eps),
   % parbolic (e=0)
   D    = RV./sqrt(Mu);
   % Mean anomaly
   M    = Q.*D + D.^3./6;
   % Mean motion
   MM   = GK.*sqrt(Mu);
else
   % e<1
   M    = EccA - AbsE.*SinE;
   % Mean motion
   MM   = GK.*sqrt(Mu./(A.^3));
end
% Time of periastron [JD]
MM   = MM./GK;
T    = JD - M./MM;
% Mean Longitude
ML   = M + (W + Omega);
ML   = (ML./(2.*pi) - floor(ML./(2.*pi))).*2.*pi;

Elem  = [Q, AbsE, T, Inc, Omega, W];
Elem2 = [A, Nu, EccA, M, ML, MM];

