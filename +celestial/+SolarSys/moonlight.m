function Illum=moonlight(Alt,HP,Elon)
% Calculate the Moon illumination in Lux on horizontal surface 
% Package: celestial.SolarSys
% Description: Calculate the Moon illumination in Lux on horizontal
%              surface as a function of the Moon altitude, horizontal
%              parallax and Elongation.
% Input  : - Vector of Altitude in radians.
%          - Vector of Horizontal Parallax in radians.
%          - Vector of Elongation in radians.
% Output : - Illumination in Lux on horiz. surface.
% See also: sunlight.m, skylight.m
% Tested : Matlab 5.2
%     By : Eran O. Ofek                    Aug 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Illum=celestial.SolarSys.moonlight(1,0.01,1)
% Reliable: 2
%------------------------------------------------------------------------------
RAD = 180./pi;

Illum = zeros(size(Alt));


for J=1:1:length(Alt(:,1))
   X = Alt(J).*RAD./90;
 
   % Handle Altitude
   if (Alt(J).*RAD>20)
      L0 = -1.95;
      L1 = 4.06;
      L2 = -4.24;
      L3 = 1.56;
      Error = 0.02;
      LI1 = L0 + L1.*X + L2.*X.*X + L3.*X.*X.*X;
   elseif (Alt(J).*RAD<=20 && Alt(J).*RAD>5)
      L0 = -2.58;
      L1 = 12.58;
      L2 = -42.58;
      L3 = 59.06;
      Error = 0.03;
      LI1 = L0 + L1.*X + L2.*X.*X + L3.*X.*X.*X;
   elseif (Alt(J).*RAD<=5 && Alt(J).*RAD>-0.8)
      L0 = -2.79;
      L1 = 24.27;
      L2 = -252.95;
      L3 = 1321.29;
      Error = 0.03;
      LI1 = L0 + L1.*X + L2.*X.*X + L3.*X.*X.*X;
   else
      % no light
      LI1 = -10;
   end

   % Handle Elongation
   f = 180 - Elon(J).*RAD;
   LI2 = -8.68e-3.*f - 2.2e-9.*f.^4;

   % Handle Parallax
   LI3 = 2.*log10(HP(J).*RAD./0.951);

   Illum(J) = 10.^(LI1+LI2+LI3);
end


