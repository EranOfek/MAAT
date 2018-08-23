function [RV,K2]=binary_rv(t,P,T,q,e,i,omega)
% Binary star radial velocity (RV)
% Package : AstroUtil.binary
% Description: Calculate binary star radial velocity as a function of
%              time.
% Input  : - Time (days).
%          - Period (days)/
%          - Time of periastron (days).
%            Default is 0.
%          - Periastron distance [q=a*(1-e)] (AU). Default is 1.
%          - Orbital eccentricity (0<=e<1). Default is 0.
%          - Orbital inclination (radians). Default is pi./2.
%          - Longitude of the periastron (radians). Default is 0.
% Output : - Radial velocity [cm/s].
%          - K2 (secondary amplitude) [AU/day].
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: t=[0:1:365].'; RV=AstroUtil.binary.binary_rv(t,365);
%          e=0.3; a=1; RV=AstroUtil.binary.binary_rv(t,365,0,a./(1-e),e,85./RAD,100);
% Reliable: 2
%--------------------------------------------------------------------------


Def.T = 0;
Def.q = 1;
Def.e = 0;
Def.i = pi./2;
Def.omega = 0;
if (nargin==2)
   T = Def.T;
   q = Def.q;
   e = Def.e;
   i = Def.i;
   omega = Def.omega;
elseif (nargin==3)
   q = Def.q;
   e = Def.e;
   i = Def.i;
   omega = Def.omega;
elseif (nargin==4)
   e = Def.e;
   i = Def.i;
   omega = Def.omega;
elseif (nargin==5)
   i = Def.i;
   omega = Def.omega;
elseif (nargin==6)
   omega = Def.omega;
elseif (nargin==7)
   % do nothing
else
    error('Illegal number of input arguments');
end

n = 2.*pi./P;
M = n.*(t-T);
[nu,R]=celestial.Kepler.kepler_elliptic(M,q,e,NaN);
a = q./(1-e);
K2 = n.*a.*sin(i)./sqrt(1-e.^2);
RV = K2.*(e.*cos(omega) + cos(nu+omega));  % AU/day
RV = RV.*constant.au./86400;
