function H=hubble_z(Z,CosmoPars)
% The Hubble parameter as a function of redshift
% Package: AstroUtil.cosmo
% Description: Compute the Hubble parameter as a function of redshift.
%              (Assuming matter dominated universe - Z<1000).
% Input  : - Vector of redshifts. [z].
%          - Cosmological parameters : [H0, \Omega_{m}, \Omega_{\Lambda}],
%            or cosmological parmeters structure, or a string containing
%            parameters source name, default is 'wmap9' (see cosmo_pars.m).
% Output : - Hubble parameter per each redshift [km/sec/Mpc].
% Reference : Lahav et al. 1991, MNRAS, 251, 128
% Tested : Matlab 5.1
%     By : Eran O. Ofek                    Jul 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: H=AstroUtil.cosmo.hubble_z(1,[70,0.3,0.7]);
% Reliable: 1
%------------------------------------------------------------------------------
C  = 29979245800;    % speed of light [cm/sec]
Pc = 3.0857e18;      % Parsec [cm]

if (nargin==1),
   CosmoPars = 'wmap9';
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguements');
end

if (ischar(CosmoPars)==0 && isstruct(CosmoPars)==0),
   % do nothing
else
   Par = AstroUtil.cosmo.cosmo_pars(CosmoPars);
   CosmoPars = [Par.H0, Par.OmegaM, Par.OmegaL, Par.OmegaRad];
end


H0     = CosmoPars(1);
OmegaM = CosmoPars(2);
OmegaL = CosmoPars(3);

% convert H0 to cm/sec/sec
H0_cgs = H0*100000/(Pc*1e6);

H = H0.*sqrt(OmegaM.*(1+Z).^3 - (OmegaM + OmegaL - 1).*(1+Z).^2 + OmegaL);
