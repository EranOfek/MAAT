function Mag=asteroid_magnitude(R,Delta,Beta,H,G)
% Calculate the magnitude of minor planets in the HG system
% Package: celestial.SolarSys
% Description: Calculate the magnitude of minor planets in the HG system.
%              Valid for phase angles (Beta) in range 0 to 120 deg.
% Input  : - MP-Sun distance in au.
%          - MP-observer distance in au.
%          - Phase angle in radians.
%          - The mean absolute visual magnitude (H).
%          - The slope parameter (G), default is 0.15.
% Output : - The minor planet visual magnitude.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Mag=celestial.SolarSys.asteroid_magnitude(3,2,0,15,0.15)
% Reliable: 2
%--------------------------------------------------------------------------
if (nargin==4),
   G = 0.15;
elseif (nargin==5),
   % do nothing
else
   error('Illigal number of input arguments');
end

Phi1 = exp(-3.33.*tan(0.5.*Beta).^0.63);
Phi2 = exp(-1.87.*tan(0.5.*Beta).^1.22);

Mag = H + 5.*log10(R.*Delta) - 2.5.*log10((1-G).*Phi1 + G.*Phi2);



