function [Nu,R,S,Vel]=kepler_parabolic(T,Q,K)
% Solve the Kepler equation for Parabolic orbit
% Package: celestial.Kepler
% Description: Solve Kepler equation (M = E - e sin E) and find the true
%              anomaly (related to the eccentric anomaly, E) and radius
%              vector for parabolic orbit (i.e., e=1). The function 
%              requires the time since periastron (t-T), the periastron 
%              distance (q) and the Gaussian constant for the system and
%              units (k; see gauss_grav_const.m).
% Input  : - Time since periastron (t-T) in units specified by chosen
%            Gaussian gravitational constant (default is days).
%          - Periastron distance (q) in units specified by chosen
%            Gaussian gravitational constant (default is au).
%          - Gaussian constant for the system (k).
%            Default is : 0.017202098950000 (appropriate for the solar
%            system; negligible mass object orbiting the Sun; time units 
%            are days and distance units are au). See gauss_grav_const.m
% Output : - True anomaly [radians].
%          - Radius vector in units of distance, for the default case units
%            are au.
%          - S (=2qS^2 = D^2) [radians].
%          - Orbital velocity [distance units per time units], for default
%            Gaussian gravitational constant this is [au/day].
% See also: kepler_elliptic.m, kepler_hyperbolic.m, kepler.m,
%           gauss_grav_const.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Nu,R,S]=celestial.Kepler.kepler_parabolic([7.0;7.8],0.1);
% Reliable: 2
%--------------------------------------------------------------------------
if (nargin==2),
   K   = 0.017202098950000;
elseif (nargin==3),
   % do nothing
else
   error('Illigal number of input arguments');
end

% parabolic motion - exact solution

N  = max([length(T), length(Q)]);
Nu = zeros(N,1);

SignT = sign(T);
T     = abs(T);

Const = 3.*K./sqrt(2);
W  = Const.*(T)./(Q.*sqrt(Q));

TanBeta  = 2./W;
Beta     = atan(TanBeta);
TanGamma = (tan(0.5.*Beta)).^(1./3);
Gamma    = atan(TanGamma);
S        = 2./(tan(2.*Gamma));

Nu = 2.*atan(S);

Nu(T==0) = 0;
Nu       = Nu.*SignT;
S        = S.*SignT;

% Radius vector
R     = Q.*(1+S.^2);


Vel = K.*sqrt(2).*sqrt(1./R);


