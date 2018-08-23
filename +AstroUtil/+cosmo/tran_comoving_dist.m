function D_M=tran_comoving_dist(Z,CosmoPars)
% Transverse comoving distance
% Package: AstroUtil.cosmo
% Description: Calculate the transverse comoving distance. Given the
%              transverse comoving distance (D_M), the comoving distance
%              between two events at the same redshift, but separated on
%              the sky by some angle \delta\theta is: D_M\delta\theta.
% Input  : - Redshift.
%          - Cosmological parameters:
%            [H0, \Omega_{m}, \Omega_{\Lambda}, \Omega_{radiation}],
%            or cosmological parmeters structure, or a string containing
%            parameters source name, default is 'wmap3' (see cosmo_pars.m).
%            Default for \Omega_{radiation} is 0.
% Output : - Comoving distance [pc].
% Reference: http://nedwww.ipac.caltech.edu/level5/Hogg/frames.html
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
% Example: D_M=AstroUtil.cosmo.tran_comoving_dist([1;2]);
%---------------------------------------------------------------------------

if (nargin==1),
   CosmoPars = 'wmap3';
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end


if (ischar(CosmoPars)==0 && isstruct(CosmoPars)==0),
   % do nothing
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

CosmoPars = [H0, OmegaM, OmegaL, OmegaRad];


H0        = CosmoPars(1);
C         = constant.c;
DH        = C./(H0.*1e5) .*1e6;               % [pc]
DC        = AstroUtil.cosmo.comoving_dist(Z,CosmoPars);       % [pc]
OmegaPars = CosmoPars(2:end);
OmegaK    = 1 - sum(OmegaPars);


if (OmegaK==0),
   D_M = DC;   
elseif (OmegaK>0),
   D_M = DH.*sinh(sqrt(OmegaK).*DC./DH)./sqrt(OmegaK);
else
   % OmegaK<0
   D_M = DH.*sin(sqrt(abs(OmegaK)).*DC./DH)./sqrt(abs(OmegaK));
end

