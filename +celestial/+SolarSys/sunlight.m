function Illum=sunlight(Alt)
% Calculate the Sun illumination in Lux on horizontal surface 
% Package: celestial.SolarSys
% Description: Calculate the Sun illumination in Lux on horizontal
%              surface as a function as its altitude in radians.
% Input  : - vector of Altitude in radians.
% Output : - Illumination in Lux on horiz. surface.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Illum=celestial.SolarSys.sunlight(-0.1)
% Reliable: 2
%------------------------------------------------------------------------------
RAD = 180./pi;

Illum = zeros(size(Alt));

for J=1:1:length(Alt(:,1)),
   X = Alt(J).*RAD./90;

   if (Alt(J).*RAD>20),
      L0 = 3.74;
      L1 = 3.97;
      L2 = -4.07;
      L3 = 1.47;
      Error = 0.02;
      LI = L0 + L1.*X + L2.*X.*X + L3.*X.*X.*X;
   elseif (Alt(J).*RAD<=20 && Alt(J).*RAD>5),
      L0 = 3.05;
      L1 = 13.28;
      L2 = -45.98;
      L3 = 64.33;
      Error = 0.02;
      LI = L0 + L1.*X + L2.*X.*X + L3.*X.*X.*X;
   elseif (Alt(J).*RAD<=5 && Alt(J).*RAD>-0.8),
      L0 = 2.88;
      L1 = 22.26;
      L2 = -207.64;
      L3 = 1034.30;
      Error = 0.02;
      LI = L0 + L1.*X + L2.*X.*X + L3.*X.*X.*X;
   elseif (Alt(J).*RAD<=-0.8 && Alt(J).*RAD>-5),
      L0 = 2.88;
      L1 = 21.81;
      L2 = -258.11;
      L3 = -858.36;
      Error = 0.02;
      LI = L0 + L1.*X + L2.*X.*X + L3.*X.*X.*X;
   elseif (Alt(J).*RAD<=-5 && Alt(J).*RAD>-12),
      L0 = 2.70;
      L1 = 12.17;
      L2 = -431.69;
      L3 = -1899.83;
      Error = 0.01;
      LI = L0 + L1.*X + L2.*X.*X + L3.*X.*X.*X;
   elseif (Alt(J).*RAD<=-12 && Alt(J).*RAD>-18),
      L0 = 13.84;
      L1 = 262.72;
      L2 = 1447.42;
      L3 = 2797.93;
      Error = 0.01;
      LI = L0 + L1.*X + L2.*X.*X + L3.*X.*X.*X;
   else
      % only starlight + airglow
      LI = -2.69897;
   end
   Illum(J) = 10.^LI;
end
