function Frac=binary_reflection_effect(R2,D,Ph,Pars)
% reflection effect from a star with unit illumination on mirror
% Package: AstroUtil.binary
% Description: Calculate the reflection effect from a star with unit
%              illumination on mirror (star) with radius R2.
% Input  : - Mirror radius [consistent length unit].
%          - Distance between stars [consistent length unit].
%          - Phase angle, observer-mirror-star, [radians].
%          - Model parameters vector:
% Output : - The fraction of reflected luminosity (from L1, by L2).
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------



SkyCov  = (R2./D).^2;           % sky coverage of L2 (in L1 system).
IllFrac = 0.5.*(1 - cos(Ph));   % Illuminated fraction

ModelType = 'Mirror';

switch ModelType
 case 'Mirror'
    % simple mirror model with Pars(1) reflectivity  
    Frac = Pars(1).*IllFrac.*SkyCov;
 otherwise
    error('Unknown reflection ModelType');
end


 
