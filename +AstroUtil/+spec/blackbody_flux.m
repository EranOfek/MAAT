function Flux=blackbody_flux(Temp,WaveRange,Radius,Dist)
% Flux of blackbody in some wavelength range
% Package: AstroUtil.spec
% Description: Calculate the flux, in a given wavelength range, of a
%              black-body, given its temperature, radius and distance.
% Input  : - Vector of black body temperature [K].
%          - Wavelength range [Min Max] in Ang.
%            Use convert.energy.m to convert from other units.
%          - Black body radius [cm], default is 1cm.
%          - Distance [pc], default is 10pc.
% Output : - Black body observed flux [erg cm^-2 s^-1].
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       Feb 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See Also: blackbody_mag.m, blackbody_bolmag.m
% Reliable: 2
%------------------------------------------------------------------------------
Pc  = constant.pc;

DefRadius = 1;
DefDist   = 10;
if (nargin==2)
   Radius = DefRadius;
   Dist   = DefDist;
elseif (nargin==3)
   Dist   = DefDist;
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

WaveRange = linspace(min(WaveRange),max(WaveRange),100).';

Nt  = length(Temp);
if (length(Radius)==1)
   Radius = Radius.*ones(Nt,1);
end
if (length(Dist)==1)
   Dist = Dist.*ones(Nt,1);
end

Flux = zeros(Nt,1).*NaN;
for It=1:1:Nt
   [Il,If] = AstroUtil.spec.black_body(Temp(It),WaveRange);
   Spec    = [WaveRange, Il.*1e-8];

   Flux(It) = trapz(Spec(:,1),Spec(:,2));
end

Flux = Flux .* 4.*pi.* Radius.^2 ./(4.*pi.*(Dist.*Pc).^2);
