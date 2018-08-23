function Dist = ad_dist(Z,CosmoPars,varargin)
% Calculate the filled beam angular diameter distance between two redshifts
% Package: AstroUtil.cosmo
% Description: Calculate the filled beam angular diameter distance
%              between two redshifts along the line of sight.
% Input  : - Two, or one columns matrix containing redshifts [z1 z2].
%            For each row, the angular diameter distance is
%            calculated between z1 and z2. If one column is given [z2],
%            then calculate the angular diameter distance from
%            z1=0 to z2.
%          - Cosmological parameters : [H0, \Omega_{m}, \Omega_{\Lambda}],
%            or cosmological parmeters structure, or a string containing
%            parameters source name, default is 'planck' (see cosmo_pars.m).
%            If a scalar is givne then assume its H_0 [km/s/Mpc].
%          - Optional Omega_{m}.
%          - Optional Omega_{\Lambda}.
% Output : - Angular diameter distance [pc].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstroUtil.cosmo.ad_dist([1 2;1.5 2],'wmap3');
% Reliable: 1
%--------------------------------------------------------------------------
Tol = 1e-4;   % integration tolerance

C  = 29979245800;    % speed of light [cm/sec]
Pc = 3.0857e18;      % Parsec [cm]

if (nargin==1),
   CosmoPars = 'planck';
end


if (ischar(CosmoPars)==0 && isstruct(CosmoPars)==0),
   % do nothing
   if (length(CosmoPars)==1),
      CosmoPars = [CosmoPars, varargin{1}, varargin{2}];
   end
else
   Par = AstroUtil.cosmo.cosmo_pars(CosmoPars);
   CosmoPars = [Par.H0, Par.OmegaM, Par.OmegaL, Par.OmegaRad];
end

H0       = CosmoPars(1);
OmegaM   = CosmoPars(2);
OmegaL   = CosmoPars(3);
if (length(CosmoPars)==3),
   OmegaRad = 0;
else
   OmegaRad = CosmoPars(4);
end


% convert H0 to cm/sec/sec
H0 = H0.*100000./(Pc.*1e6);

N = size(Z,1);
if (size(Z,2)==1),
   Z1 = zeros(N,1);
   Z2 = Z;
else
   Z1 = Z(:,1);
   Z2 = Z(:,2);
end

% calculate general geometry
OmegaTot  = OmegaM + OmegaL;

R0 = C./H0;  % cm.

OmegaK = 1 - OmegaM - OmegaL - OmegaRad;
Ez = @(Z) 1./sqrt(OmegaRad.*(1+Z).^4 + OmegaM.*(1+Z).^3 + OmegaK.*(1+Z).^2 + OmegaL);


IntVal = zeros(N,1);
for I=1:1:N,
   % integrate inv_e_z from z1 to z2
   %IntVal(I) = quad(inv_e_z',Z1(I),Z2(I),Tol,[],CosmoPars(2:end));
   IntVal(I) = integral(Ez,Z1(I),Z2(I),'AbsTol',Tol); %,Tol,[],CosmoPars(2:end));
end

if (OmegaTot>1),
   % k=+1   close
   Chi12 = sqrt(abs(OmegaTot - 1)).*IntVal;
   Dist  = R0.*sin(Chi12)./((1+Z2).*sqrt(OmegaTot-1));
elseif (OmegaTot<1),
   % k=-1   open
   Chi12 = sqrt(abs(OmegaTot - 1)).*IntVal;
   Dist  = R0.*sinh(Chi12)./((1+Z2).*sqrt(1-OmegaTot));
else
   % k=0    flat
   Chi12 = IntVal;
   Dist  = R0.*Chi12./(1+Z2);
end


% convert distance in cm to parsecs
Dist = Dist./Pc;





