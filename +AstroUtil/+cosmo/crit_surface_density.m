function SigmaCrit=crit_surface_density(Z,CosmoPars,DistType)
% The critical surface density for gravitational lensing
% Package: AstroUtil.cosmo
% Description: Calculates the critical surface density for gravitational
%              lensing, given the cosmology and redshifts of the lens
%              and source.
% Input  : - Redshifts vectors [[z_{l}], [z_{s}]], or distances in parsecs.
%          - Cosmological parameters : [H0, \Omega_{m}, \Omega_{\Lambda}],
%            or cosmological parmeters structure, or a string containing
%            parameters source name, default is 'wmap3' (see cosmo_pars.m).
%          - Distances type:
%            'norm' : normal angular diam. dist. - default.
%            'dr'   : Dyer-Roeder ang. diam. distance (empty beam case).
%            'pc'   : flat space (nearby sources) - In this case the input
%                     is distances in parsecs (instead of redshift).
% Output : - Critical surface density for gravitational lensing [gr/cm^2].
% Tested : Matlab 5.1
%     By : Eran O. Ofek                    Jul 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: SigmaCrit=AstroUtil.cosmo.crit_surface_density([1 2],'wmap3');
% Reliable: 2
%--------------------------------------------------------------------------
G = 6.672e-8;
C = 29979245800;
Pc = 3.0857e18;

if (nargin==1),
   CosmoPars = 'wmap3';
   DistType  = 'norm';
elseif (nargin==2)
   DistType  = 'norm';
elseif (nargin==3),
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


Zl = Z(:,1);
Zs = Z(:,2);

switch DistType
 case 'pc'
    % do nothing
 otherwise
    H0     = CosmoPars(1);
    OmegaM = CosmoPars(2);
    OmegaL = CosmoPars(3);
end

switch DistType
 case 'norm'
    %Dz = inline('ad_dist(Z,H0,OmegaM,OmegaL)','Z','H0','OmegaM','OmegaL');
    Dz = @(Z,H0,OmegaM,OmegaL) ad_dist(Z,H0,OmegaM,OmegaL);
 case 'dr'
    %Dz = inline('ad_dr_dist(Z,H0,OmegaM,OmegaL)','Z','H0','OmegaM','OmegaL');
    Dz = @(Z,H0,OmegaM,OmegaL) ad_dr_dist(Z,H0,OmegaM,OmegaL);
 case 'pc'
    % do nothing
 otherwise
    error('Unknown DistType option');
end

switch DistType
 case 'pc'
    Ds  = Zs.*Pc;
    Dl  = Zl.*Pc;
    Dls = (Zs - Zl).*Pc;
 otherwise
    Ds  = Dz(Zs,H0,OmegaM,OmegaL).*Pc;
    Dl  = Dz(Zl,H0,OmegaM,OmegaL).*Pc;
    Dls = Dz([Zl, Zs],H0,OmegaM,OmegaL).*Pc;
end

SigmaCrit = (C.^2./(4.*pi.*G)).*Ds./(Dl.*Dls);
