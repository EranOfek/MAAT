function T=lookback_time(Zvec,CosmoPars)
% Compute the cosmological lookback time
% Package: AstroUtil.cosmo
% Description: Compute the cosmological lookback time, between two events
%              in redshift z1 and z2, and given the cosmology.
%              (Assuming matter dominated universe - Z<1000).
% Input  : - Matrix of redshifts. [[z1], [z2]], if only one column
%            is given then takes z1=0.
%          - Cosmological parameters : [H0, \Omega_{m}, \Omega_{\Lambda}],
%            or cosmological parmeters structure, or a string containing
%            parameters source name, default is 'wmap3' (see cosmo_pars.m).
%            Ignores \Omega_{\radiation}.
% Output : - Lookback time [seconds].
% Reference : Lahav et al. 1991, MNRAS, 251, 128
% Tested : Matlab 5.1
%     By : Eran O. Ofek                    Jul 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: T=AstroUtil.cosmo.lookback_time([0 1; 1 2; 2 3],[70 0.1 0.9])./(3.1e7.*1e9)
% Reliable: 1
%------------------------------------------------------------------------------
C  = 29979245800;    % speed of light [cm/sec]
Pc = 3.0857e18;      % Parsec [cm]

if (nargin==1),
   CosmoPars = [70, 0.3, 0.7];
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

SizeZ = size(Zvec);
if (SizeZ(2)==1),
   Z1 = zeros(SizeZ(1),1);
   Z2 = Zvec;
else
   Z1 = Zvec(:,1);
   Z2 = Zvec(:,2);
end


A_high = 1./(1+Z1);
A_low  = 1./(1+Z2);

% lookback time integrand (still need to divide by H0)
%lookback_fun = inline('sqrt(A./(OmegaL.*A.^3 + (1-OmegaM-OmegaL).*A + OmegaM))','A','OmegaM','OmegaL');
lookback_fun = @(A,OmegaM,OmegaL) sqrt(A./(OmegaL.*A.^3 + (1-OmegaM-OmegaL).*A + OmegaM));


T = zeros(SizeZ(1),1);
for I=1:1:SizeZ(1),
   T(I) = quadl(lookback_fun,A_low(I),A_high(I),[],[],OmegaM,OmegaL);
   T(I) = T(I)./H0_cgs;
end

