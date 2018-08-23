function [Nu,R,H,Vel]=kepler_hyperbolic(T,Q,E,K,Tol)
% SOlve Kepler equation for hyperpolic orbit
% Package: celestial.Kepler
% Description: Solve Kepler equation (M = e sinh H - H) and find the true
%              anomaly and radius vector for hyperbolic orbit (i.e., e>1).
%              The function requires the time since periastron (t-T), the
%              periastron distance (q) and the Gaussian constant for the
%              system and units (k; see gauss_grav_const.m).
% Input  : - Time since periastron (t-T) in units specified by chosen
%            Gaussian gravitational constant (default is days).
%          - Periastron distance (q) in units specified by chosen
%            Gaussian gravitational constant (default is au).
%          - Orbital eccentricity.
%          - Gaussian constant for the system (k).
%            Default is : 0.017202098950000 (appropriate for the solar
%            system; negligible mass object orbiting the Sun; time units
%            are days and distance units are au). See gauss_grav_const.m
%          - Tolerance, default is 1e-8.
% Output : - True anomaly [radians].
%          - Radius vector in units of distance, for the default case units
%            are au.
%          - Hyperbolic eccentric anomaly, H [radians].
%          - Orbital velocity [distance units per time units], for default
%            Gaussian gravitational constant this is [au/day].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Nu,R,H,Vel]=celestial.Kepler.kepler_hyperbolic(1,1,1.1);
% Reliable: 2
%--------------------------------------------------------------------------
DefTol = 1e-8;
if (nargin==3),
   K   = 0.017202098950000;
   Tol = DefTol;
elseif (nargin==4),
   Tol = DefTol;
elseif (nargin==5),
   % do nothing
else
   error('Illigal number of input arguments');
end


MaxLength = max([length(T), length(Q), length(E)]);
if (length(T)==1),
  T = T.*ones(MaxLength,1);
end
if (length(Q)==1),
    Q = Q.*ones(MaxLength,1);
end
if (length(E)==1),
    E = E.*ones(MaxLength,1);
end



% hyperbolic motion

SignT = sign(T);
T     = abs(T);
% q = a(1 - e), a is negative
A = Q./(1 - E);
MM = K./sqrt((-A).^3);  % mean motion
M  = MM.*T;             % mean anomaly

H0 = M;
H1 = H0 + (M - E.*sinh(H0) + H0)./(-1 + E.*cosh(H0));

Ns = length(H1);
IndUnsolved = (1:1:Ns).';                                                       

I = 0;
%while max(abs(H1(IndUnsolved)-H0(IndUnsolved)))>Tol,
while (any(abs(H1-H0)>Tol)),
       
   I = I + 1;
   H0 = H1;
   H1 = H0 + (M - E.*sinh(H0) + H0)./(-1 + E.*cosh(H0));
%    H0(IndUnsolved) = H1(IndUnsolved);
%    H1(IndUnsolved) = H0(IndUnsolved) + ...
%                     (M(IndUnsolved) - E(IndUnsolved).*sinh(H0(IndUnsolved)) + H0(IndUnsolved))./(-1 + E(IndUnsolved).*cosh(H0(IndUnsolved)));
%    IndUnsolved = find(abs(H1(IndUnsolved)-H0(IndUnsolved))>Tol);  
end
H = H1;

CosNu = (cosh(H) - E)./(1 - E.*cosh(H));
Nu    = acos(CosNu);

H     = H.*SignT;
Nu    = Nu.*SignT;


%R = Q.*(1+E)./(1+E.*cos(Nu));
R = Q.*(1+E)./(1+E.*CosNu);


A = Q./(1-E);
Vel = K.*sqrt(2).*sqrt(1./R - 1./(2.*A));


