function DM=disp_measure(Z,Par)
% Calculate the dispersion measure along a cosmological line of sight
% Package: AstroUtil.cosmo
% Description: Calculate the dispersion measure along a cosmological line
%              of sight as a function of redshift. Valid below redshift 3
%              (He re-ionization).
% Input  : - Vector of redshifts.
%          - An optional structure of cosmological parameters (see
%            AstroUtil.cosmo.cosmo_pars). Alternatively a string of
%            parameters names. Default is 'planck'.
% Output : - Vector of dispersion measures [cm^-3 pc].
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: Zheng, Z. et al. 2014 ApJ 797, 71
% Example: DM=AstroUtil.cosmo.disp_measure([0.1;1])
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<2)
    Par = 'planck';
end
if (isstruct(Par))
    % do nothing
elseif (ischar(Par))
    Par = AstroUtil.cosmo.cosmo_pars(Par);
else
    error('Illegal cosmological parameters input');
end

OmegaM  = Par.OmegaM;
OmegaL  = Par.OmegaL;

Fe = 0.88; % below z~3
Fun = @(z) (1+z).*Fe./sqrt(OmegaM.*(1+z).^3 + OmegaL);

Nz = numel(Z);
DM = zeros(Nz,1);
for Iz=1:1:Nz
    DM(Iz) = 1060.*(Par.OmegaB_h2./0.022).*(Par.H0./70).^-1.*integral(Fun,0,Z(Iz));
end
