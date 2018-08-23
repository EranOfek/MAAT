function [LorentzFactorA,LorentzFactorB]=lorentz_from_flux(Flux,DeltaT,Z,Dist,E_min,E_max,Alpha)
% Lower limit on Lorentz factor GRB
% Package: AstroUtil.GRB
% Description: Given a GRB flux, estimate a lower limit on the Lorentz
%              factor og a GRB, assuming a power law spectrum.
% Input  : - Flux [erg/cm^2/s]
%          - DeltaT [s]
%          - Redshift
%          - Luminosity Dist [pc] (if empty use redshift).
%          - E_min [keV]
%          - E_max [keV]
%          - Alpha (spectral index).
% Output : - Lower limit on lorentz factor assuming photon can
%            escape without annihilating other photons.
%            (case A in Lithwick & Sari).
%          - Lower limit on lorentz factor assuming e+ e- pairs
%            produced by photin annihilation are optically thin
%            (case B in Lithwick & Sari).
% Reference : Lithwick & Sari 2001 ApJ 555, 540L
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Aug 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

if (isempty(Dist)==1),
   Dist = lum_dist(Z);
end

SigmaT   = constant.sigmaT;
C        = constant.c;
Me       = constant.me;
Pc       = constant.pc;

Dist = Dist.*Pc;

E_max=convert.energy('eV','erg',E_max.*1e3);
E_min=convert.energy('eV','erg',E_min.*1e3);
Int = (E_max.^(-Alpha+2) - E_min.^(-Alpha+2))./(-Alpha+2);
F = Flux./Int;

% all the photons in Epeak...
%F = Flux./convert_energy((1920).*1e3,'eV','erg') .* convert_energy(1920e3,'eV','erg').^Alpha

TauHat   = (11./180).*SigmaT.*Dist.^2.*(Me.*C.^2).^(-Alpha+1) .* F./(C.^2 .* DeltaT.*(Alpha-1));

LorentzFactorA = TauHat.^(1./(2.*Alpha+2)) .* (E_max./(Me.*C.^2)).^((Alpha-1)./(2.*Alpha+2)) .* (1+Z).^((Alpha-1)./(Alpha+1));

LorentzFactorB = TauHat.^(1./(Alpha+3)) .* (1+Z).^((Alpha-1)./(Alpha+3));



