function Ez=e_z(Z,CosmoPars)
% E(z) cosmological function
% Package: AstroUtil.cosmo
% Description: Calculate E(z) cosmological function, which is
%              proportional to the time derivative of the logarithm
%              of the scale factor.
% Input  : - Redshift.
%          - Cosmological parameters vector:
%            [Omega_M Omega_Lambda Omega_Radiation],
%            default for Omega_radiation is 0.
% Output : - E(z).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Ez=AstroUtil.cosmo.e_z([1;2],[0.3 0.7]);
% Reliable: 2
%--------------------------------------------------------------------------

if (size(CosmoPars,2)==2),
   CosmoPars = [CosmoPars, 0];
end

OmegaM = CosmoPars(1);
OmegaL = CosmoPars(2);
OmegaR = CosmoPars(3);
OmegaK = 1 - OmegaM - OmegaL - OmegaR;

Ez = sqrt(OmegaR.*(1+Z).^4 + OmegaM.*(1+Z).^3 + OmegaK.*(1+Z).^2 + OmegaL);

