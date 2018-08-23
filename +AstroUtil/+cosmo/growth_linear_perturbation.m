function [D,Z]=growth_linear_perturbation(Z,OmegaPar)
% alculate the growth function of linear perurbations
% Package: AstroUtil.cosmo
% Description: Calculate the growth function of linear perurbations in
%              various cosmological models.
% Input  : - Vector of redshifts.
%          - Omega parameters, [\Omega_{m}, \Omega_{\Lambda}].
% Output : - Growth speed vector for each Z.
%            normalized to 1 at z=0.
%          - Vector of input redshifts.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: Li & Ostriker, 2000 astro-ph/0010432
%            Peebles 1980, The large scale structure of the universe.
% Example: [D]=AstroUtil.cosmo.growth_linear_perturbation(1,[0.3 0.7]);
%--------------------------------------------------------------------------

OmegaM = OmegaPar(1);
OmegaL = OmegaPar(2);
OmegaK = 1 - OmegaM - OmegaL;

if (OmegaM==1 && OmegaL==0),
   % Einstein-de Sitter model
   D = 1./(1+Z);
elseif (OmegaM>0 && OmegaM<1 && OmegaL==0),
    % Open model
    %F1 = inline('1+3./x + 3.*sqrt(1+x).*log(sqrt(1+x) - sqrt(x))./(x.^1.5)','x');
    F1 = @(x) 1+3./x + 3.*sqrt(1+x).*log(sqrt(1+x) - sqrt(x))./(x.^1.5);
    D = F1((1./OmegaM - 1)./(1+Z))./F1(1./OmegaM - 1);
elseif (OmegaM>0 && OmegaM<1 && OmegaL==(1-OmegaM) && OmegaL<1),
   % \Lambda-model (with errors <1%)
   %F2 = inline('0.358.*(1-(1+0.23.*x.^2.2).^-1.26).^(1./2.2)','x');
   F2 = @(x) 0.358.*(1-(1+0.23.*x.^2.2).^-1.26).^(1./2.2);
   Arg1 = ((2.*OmegaL/OmegaM).^(1./3))./(1+Z);
   Arg2 = (2.*OmegaL/OmegaM).^(1./3);
   D = F2(Arg1)./F2(Arg2);
else
   error('Unkown model');
end

