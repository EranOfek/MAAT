function [I_rho,A_rho]=fresnel_occultation_ps(Rho,r)
% Monchromatic diffraction pattern for finite source
% Package: AstroUtil.Occultation
% Description: Calculate the intensity pattern due to
%              diffractive Fresnel occultation by a small circular object.
%              Assuming a mono-chromatic observation and a point source.
% Input  : - Radius (scalar) of the occulting object in units of the
%            Fresnel scale, where the Fresnel scale is \sqrt(d\lambda/2),
%            where \lambda is the wavelength in which the observations
%            are conducted and d is the distance from the observer to
%            the occulting object.
%          - Vector of distances of the center of the occulting object
%            from the line of sight to the source, in units of the
%            Fresnel scale.
% Output : - Intensity as a function of the distance from the line of
%            sight.
%          - Complex intensity and phase as a function of the distance
%            from the line of sight.
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Jun 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: fresnel_occultation_mono.m
% Example: [I_rho,A_rho]=AstroUtil.Occultation.fresnel_occultation_ps(0.9,[0:0.1:20]');
% Reliable: 2
%--------------------------------------------------------------------------

% calculate diffractive occultation

%FunInt = inline('exp(0.5.*i.*pi.*R.^2).*besselj(0,pi.*r.*R).*R','R','r');
FunInt = @(R,r) exp(0.5.*1i.*pi.*R.^2).*besselj(0,pi.*r.*R).*R;

Nr = length(r);
   
%Is = find(r<10);
%Il = find(r>=10);
   
%Fr     = zeros(Nr,1);
%Fr(Is) = quadv(FunInt,0,Rho,[],[],r(Is));
%Fr(Il) = 0;
   
Fr = quadv(FunInt,0,Rho,[],[],r);
   
%Fr = zeros(Nr,1);
%for Ir=1:1:Nr,
%   Fr(Ir) = quad(FunInt,0,Rho,[],[],r(Ir));
%end
   
A_rho = 1 + 1i.*pi.*exp(0.5.*1i.*pi.*r.^2).*Fr;
   
I_rho = A_rho.*conj(A_rho);

%else
%   % Assume geometric occultation
%
%   I_rho = zeros(size(r));
%
%   Ilr = find(r>Rho);
%   Isr = find(r<=Rho);
%
%   I_rho(Isr) = 0;
%   I_rho(Ilr) = 1;
%   A_rho = I_rho;
%end


