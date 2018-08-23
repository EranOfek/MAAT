function [Z]=inv_comoving_volume(Volume,CosmoPars)
% Convert cosmological volume to redshift
% Package: AStroUtil.cosmo
% Description: Use the cosmological volume to calculate the corresponding
%              redshift.
% Input  : - Vector of volumes [pc^3].
%          - Cosmological parameters : [H0, \Omega_{m}, \Omega_{\Lambda}],
%            or cosmological parmeters structure, or a string containing
%            parameters source name, default is 'wmap9' (see cosmo_pars.m).
% Output : - Redshift.
% Reference : http://nedwww.ipac.caltech.edu/level5/Hogg/frames.html
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Dec 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [~,V]=AstroUtil.cosmo.comoving_volume(2); [Z]=AstroUtil.cosmo.inv_comoving_volume(V)
% Reliable: 2
%---------------------------------------------------------------------------
Def.Type = 'LD';
Def.CosmoPars = 'wmap9';

if (nargin==1)
   Type      = Def.Type;
   CosmoPars = Def.CosmoPars;
elseif (nargin==2)
   CosmoPars = Def.CosmoPars;
elseif (nargin==5)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (ischar(CosmoPars)==0 && isstruct(CosmoPars)==0)
   % do nothing
else
   Par = AstroUtil.cosmo.cosmo_pars(CosmoPars);
   CosmoPars = [Par.H0, Par.OmegaM, Par.OmegaL, Par.OmegaRad];
end
H0       = CosmoPars(1);
OmegaM   = CosmoPars(2);
OmegaL   = CosmoPars(3);
if (length(CosmoPars)==3)
   OmegaRad = 0;
else
   OmegaRad = CosmoPars(4);
end

Ndm     = numel(Volume);
Z       = zeros(Ndm,1);
for Idm=1:1:Ndm
   Z(Idm) = Util.find.fun_binsearch(@comoving_volume_handle,Volume(Idm),[1e-15 100],1e-3,CosmoPars);
end

function Vol=comoving_volume_handle(z,CosmoPars)
[~,Vol]=AstroUtil.cosmo.comoving_volume(z,CosmoPars);


