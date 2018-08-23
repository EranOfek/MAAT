function Diff=cdt_dz(Z,CosmoPars)
% Calculate the differential cdt/dz in the FLRW geometry.
% Package: AstroUtil.cosmo
% Description: Calculate the differential cdt/dz in the FLRW geometry.
% Input  : - Vector of redshifts.
%          - Cosmological parameters: [H0, \Omega_{m}, \Omega_{\Lambda}],
%            or cosmological parmeters structure, or a string containing
%            parameters source name, default is 'wmap3' (see cosmo_pars.m).
% Output : - The quantity cdt/dz in the FLRW geometry [cm].
% Reference : Fukugita et al. 1992, ApJ 393, 3
% Tested : Matlab 5.1
%     By : Eran O. Ofek             Jul 2001
%    URL : http://weizmann.ac.il/home/eofek/
% Example: Diff=AstroUtil.cosmo.cdt_dz([1;2]);
% Reliable: 1
%------------------------------------------------------------------------------
C = 29979245800;  % speed of light [cm/sec]
Pc = 3.0857e18;   % parsec [cm]

if (nargin==1),
   CosmoPars = [70, 0.3, 0.7];
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end


if (ischar(CosmoPars)==0 && isstruct(CosmoPars)==0),
   % do nothing
else
   Par = AStroUtil.cosmo.cosmo_pars(CosmoPars);
   CosmoPars = [Par.H0, Par.OmegaM, Par.OmegaL, Par.OmegaRad];
end
if (length(CosmoPars)==3),
   CosmoPars(4) = 0;
end

H0     = CosmoPars(1);
OmegaM = CosmoPars(2);
OmegaL = CosmoPars(3);

% convert H0 to cgs:
H0 = H0.*100000./(1e6.*Pc);

% the Hubble distance
R0 = C/H0;

X    = (1 + Z);
Diff = R0./(X.*sqrt(OmegaM.*X.^3 + (1 - OmegaM - OmegaL).*X.^2 + OmegaL));





