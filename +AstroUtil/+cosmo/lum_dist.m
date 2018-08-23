function [DL,DM]=lum_dist(Z,CosmoPars)
% Luminosity distance
% Package: AstroUtil.cosmo
% Description: Compute luminosity distance from redshift and cosmological
%              parameters. Given the object spectra, calculate also the
%              K-correction.
% Input  : - Vector of redshifts.
%          - Cosmological parameters : [H0, \Omega_{m}, \Omega_{\Lambda}],
%            or cosmological parmeters structure, or a string containing
%            parameters source name, default is 'wmap3' (see cosmo_pars.m).
% Output : - Luminosity distance [parsec].
%          - Distance modulus [mag].
% Reference : Perlmutter et al. 1997 ApJ, 483, 565
%             Oke & Sandage 1968 ApJ, 154, 21
%             Peterson, B.M., 1997, AGN, p.165
% Tested : Matlab 5.1
%     By : Eran O. Ofek                    Jul 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DL,DM]=AstroUtil.cosmo.lum_dist([0.1;0.2])
%          [DL,DM]=AstroUtil.cosmo.lum_dist([0.1;0.2],[71 0.3 0.7])
% Reliable: 2
%------------------------------------------------------------------------------
C  = 29979245800;    % speed of light [cm/sec]
Pc = 3.0857e18;      % Parsec [cm]

if (nargin==1),
   CosmoPars = 'wmap3';
   Lambda = NaN;
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

% convert H0 to cm/sec/sec
H0 = H0.*100000./(Pc.*1e6);



if ((OmegaM+OmegaL)>1),
   %Lx = inline('sin(x)','x');
   Lx = @(x) sin(x);
   K  = 1 - OmegaM - OmegaL;
elseif ((OmegaM+OmegaL)<1),
   %Lx = inline('sinh(x)','x');
   Lx = @(x) sinh(x);
   K  = 1 - OmegaM - OmegaL;
else
   % OmegaM + OmegaL == 1
   %Lx = inline('x','x');
   Lx = @(x) x;
   K  = 1;
end

N  = length(Z);
DL = zeros(size(Z));
for I=1:1:N,
   ZI = Z(I);

   %Int = inline('((1+ZI).^2.*(1+OmegaM.*ZI) - ZI.*(2+ZI).*OmegaL).^(-1./2)','ZI','OmegaM','OmegaL');
   Int = @(ZI) ((1+ZI).^2.*(1+OmegaM.*ZI) - ZI.*(2+ZI).*OmegaL).^(-1./2);
   
   %DL_Int = quad(Int,0,ZI,[],[],OmegaM,OmegaL);
   DL_Int = integral(Int,0,ZI); 

   DL(I) = (C.*(1+ZI)./(H0.*sqrt(abs(K)))).*Lx(sqrt(abs(K)).*DL_Int);

end

% luminosity distance in Parsecs:
DL = DL.'./Pc;

% distance modulus:
if (nargout>1),
   DM = 5.*log10(DL./10);
   DM = DM.';
end
DL = DL.';

