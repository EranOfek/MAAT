function [Nu,R,E,Vel,M]=kepler_elliptic(T,Q,Ecc,K,Tol)
% Solve Kepler equation for elliptic orbit
% Package: celestial.Kepler
% Description: Solve Kepler equation (M = E - e sin E) and find the true
%              anomaly and radius vector for elliptical orbit (i.e., 0<=e<1).
%              The function requires the time since periastron (t-T), the
%              periastron distance (q) and the Gaussian constant for the
%              system and units (k; see gauss_grav_const.m).
% Input  : - Time since periastron (t-T) in units specified by chosen
%            Gaussian gravitational constant (default is days).
%            If the Gaussian gravitational constant is NaN, then assumes
%            this is the mean anomaly.
%          - Periastron distance (q) in units specified by chosen
%            Gaussian gravitational constant (default is au).
%            The periastron distance is related to the semi major
%            major axis (a), through q=a*(1-e).
%          - Orbital eccentricity.
%          - [K,M] Gaussian constant for the system (k) and Mass [SolarMass].
%            Default is : 
%            k=0.017202098950000 , M=1 [SolarMass];
%            (appropriate for the solar system; 
%            negligible mass object orbiting the Sun; time units
%            are days and distance units are au). See gauss_grav_const.m
%            If NaN, then use the first argument as the mean anomaly, 
%            instead of time since periastron.
%          - Tolerance, default is 1e-8.
% Output : - True anomaly [radians].
%          - Radius vector in units of distance, for the default case units
%            are au.
%          - Eccentric anomaly [radians].
%          - Orbital velocity [distance units per time units], for default
%            Gaussian gravitational constant this is [au/day].
%          - Velocity is not defined if Gaussian grav. constant is NaN.
%          - The mean anomaly [radians].
% Notes  : The mean anomaly (M) is related to the time since periastron
%           and q through: M = n*(t-T), where the mean motion (n) is
%           given by: k*a^(3/2), where a=q/(1-e) and k is the Gaussian
%           gravitational constant.
% See also: kepler_parabolic.m, kepler_hyperbolic.m, gauss_grav_const.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Nu,R,E]=celestial.Kepler.kepler_elliptic(10,0.5,0.96);
%          [Nu,R,E]=celestial.Kepler.kepler_elliptic(1, 0.5,0.96,NaN);
% Reliable: 2
%--------------------------------------------------------------------------
DefTol = 1e-8;
if (nargin==3)
   K   = 0.017202098950000;
   Mass=1;
   Tol=DefTol;
elseif (nargin==4)
    if (numel(K) == 1)
        Mass = 1; 
    else
        Mass=K(2);
        K=K(1);
    end
    Tol=DefTol;

elseif (nargin==5)
    if (numel(K) == 1)
        Mass = 1; 
    else
        Mass=K(2);
        K=K(1);
    end
    % do nothing
   Tol = DefTol;
elseif (nargin==4)
   Tol = DefTol;
elseif (nargin==5)
   % do nothing
else
   error('Illigal number of input arguments');
end
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

MaxLength = max([length(T), length(Q), length(Ecc)]);
if (length(T)==1)
   T = T.*ones(MaxLength,1);
end
if (length(Q)==1)
   Q = Q.*ones(MaxLength,1);
end
if (length(Ecc)==1)
   Ecc = Ecc.*ones(MaxLength,1);
end


% elliptic motion
if (~isnan(K))
   A      = Q./(1-Ecc);
   %PeriodEran = 2.*pi.*(A.^1.5)./K;   % in time units - true for M= 1 solar mass
   %Periodold = 2.*pi.*(A.^1.5)./K./sqrt(Mass);
   %Period = 2*pi*sqrt(((A.*(constant.au)).^3)./((constant.G).*((Mass)).*(constant.SunM)))./86400; %Period [days] 
   Period = 2.*pi.*sqrt( (A.*(constant.au)).^3 ./ (constant.G.*Mass.*constant.SunM) )./86400;  % Period [days]
   
   N      = 2.*pi./Period;        % mean motion (n)
   M      = N.*T;                 % Mean anomaly n*(t-T)
   M      = mod(M,2.*pi);
else
   % use first argument as mean anomaly
   M      = T;
end


E0 = pi.*ones(size(M));
E1 = E0 + (M + Ecc.*sin(E0) - E0)./(1 - Ecc.*cos(E0));

Ns = length(E1);
IndUnsolved = (1:1:Ns).';

I = 0;
%while max(abs( E1(IndUnsolved)-E0(IndUnsolved) ))>Tol,
while (any(abs(E1-E0)>Tol))
   I = I + 1;
%   [I, length(IndUnsolved)]

%    E0(IndUnsolved) = E1(IndUnsolved);
%    E1(IndUnsolved) = E0(IndUnsolved) + ...
%                     (M(IndUnsolved) + Ecc(IndUnsolved).*sin(E0(IndUnsolved)) - E0(IndUnsolved))./(1 - Ecc(IndUnsolved).*cos(E0(IndUnsolved)));

    E0 = E1;
    E1 = E0 + (M + Ecc.*sin(E0) - E0)./(1 - Ecc.*cos(E0));

     
   % The next line slow down the code:
   %IndUnsolved = find(abs(E1(IndUnsolved)-E0(IndUnsolved))>Tol);
end
E = E1;

% E = angle_in2pi(E);
% Di = M - (E - Ecc.*sin(E));
% max(abs(sin(Di)))
% 


TanVH = sqrt((1+Ecc)./(1-Ecc)).*tan(0.5.*E);
Nu    = 2.*atan(TanVH);

R = Q.*(1+Ecc)./(1+Ecc.*cos(Nu));


A   = Q./(1-Ecc);
Vel = K.*sqrt(2).*sqrt(1./R - 1./(2.*A));

%dNudt = N.*sqrt(1-Ecc.^2)./((1-Ecc.*cos(E)).^2);

