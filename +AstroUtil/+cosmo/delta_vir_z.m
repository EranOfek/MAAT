function [DeltaVirZ]=delta_vir_z(Z,OmegaPar)
% Calculate the virial overdensity \Delta_{vir}
% Package: AStroUtil.cosmo
% Description: Calculate the virial overdensity \Delta_{vir}, as a
%              function of redshift z, and cosmological paramaeters.
% Input  : - Redshift vector.
%          - Parameters at zero redshift.
%            [\Omega_{m}, \Omega_{\Lambda}]
% Output : - Vector of \Delta_{vir}(Z) corresponding to redshift vector.
% Reference : Bryan & Norman 1998 ApJ 495, 80
%             Bullock et al. 2000 - astro-ph/9908159
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DeltaVirZ]=AstroUtil.cosmo.delta_vir_z(0,[0.3 0.7]);
% Reliable: 2
%------------------------------------------------------------------------------

OmegaM = OmegaPar(1);
OmegaL = OmegaPar(2);

[OmegaMz]=AstroUtil.cosmo.omega_z(Z,OmegaPar);
X = OmegaMz - 1;

if (OmegaL==0),
   DeltaVirZ = (18.*pi.*pi + 60.*X - 32.*X.*X)./OmegaMz;
else
   DeltaVirZ = (18.*pi.*pi + 82.*X - 39.*X.*X)./OmegaMz;
end
