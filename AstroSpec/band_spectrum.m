function [N_EE,N_E]=band_spectrum(E,E0,Alpha,Beta);
%------------------------------------------------------------------------------
% band_spectrum function                                             AstroSpec
% Description: Calculate an un-normalized Band-spectrum (Band et al. 1993).
% Input  : - Vector of Energy [keV].
%          - E0 [keV].
%          - Alpha, default is -1.
%          - Beta, defaukt is -2.
% Output : - Number of photons (unnormalized) per cm^2/s/keV
%            To normalize spectrum, multiply by a normalization factor.
%          - Energy [keV] (unnormalized) per cm^2/s/keV
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       Feb 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

DefAlpha  = -1;
DefBeta   = -2;
if (nargin==2),
   Alpha   = DefAlpha;
   Beta    = DefBeta;
elseif (nargin==3),
   Beta    = DefBeta;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end

N_EE = zeros(size(E));
Il   = find(E>((Alpha-Beta).*E0));
Is   = find(E<=((Alpha-Beta).*E0));

N_EE(Il) = ((Alpha-Beta).*E0./100).^(Alpha-Beta) .* exp(Beta-Alpha).*(E(Il)./100).^Beta;
N_EE(Is) = (E(Is)./100).^Alpha .* exp(-E(Is)./E0);


N_E = N_EE.*E;

