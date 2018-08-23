function [Tran,Tau,CrossSec]=xray_abs(Energy,ColDen)
%--------------------------------------------------------------------------
% xray_abs function                                              AstroSpec
% Description: Given the neutral Hydrogen column density, calculate the
%              bound-free attenuation of X-rays as a function of wavelength
%              in the ISM.
%              The program assumes abundences from Ebihara (1982).
%              Absorption is due to neutral species only.
%              Adopted from Zombeck (1990).
% Input  : - Energy [keV] in range 0.03 to 10.0 keV.
%          - Column density [cm^-2].
% Output : - The fraction of transmitted X-rays (1-no attanuation,
%                                                0-full attanuation)
%          - Optical depth.
%          - Cross section.
% Tested : Matalb 7.3
%     By : Eran O. Ofek                    Nov 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: Zombeck (1990)
%            http://ads.harvard.edu/cgi-bin/bbrowse?book=hsaa&chap=6&page=0202
%            From: Morrison & McCammon (1983; ApJ 270, 119)
% Example: [Tran,Tau,CrossSec]=xray_abs([0.03:0.01:10]',1e20);
% Reliable: 2
%--------------------------------------------------------------------------
Eps = 10.*eps;

%  Energy(keV)   C0     C1      C2
D=[0.030        17.3   608.1  -2150
   0.100-Eps    17.3   608.1  -2150
   0.100        34.6   267.9   -476.1
   0.284-Eps    34.6   267.9   -476.1
   0.284        78.1    18.8      4.3
   0.400-Eps    78.1    18.8      4.3
   0.400        71.4    66.8    -51.4
   0.532-Eps    71.4    66.8    -51.4
   0.532        95.5   145.8    -61.1
   0.707-Eps    95.5   145.8    -61.1
   0.707       308.9  -380.6    294.0
   0.867-Eps   308.9  -380.6    294.0
   0.867       120.6   169.3    -47.7
   1.303-Eps   120.6   169.3    -47.7
   1.303       141.3   146.8    -31.5
   1.840-Eps   141.3   146.8    -31.5
   1.840       202.7   104.7    -17.0
   2.471-Eps   202.7   104.7    -17.0
   2.471       342.7    18.7      0.0
   3.210-Eps   342.7    18.7      0.0
   3.210       352.2    18.7      0.0
   4.038-Eps   352.2    18.7      0.0
   4.038       433.9    -2.4      0.75
   7.111-Eps   433.9    -2.4      0.75
   7.111       629.0    30.9      0.0
   8.331-Eps   629.0    30.9      0.0
   8.331       701.2    25.2      0.0
   10.00-Eps   701.2    25.2      0.0];

%Ne = length(Energy);

C0 = interp1(D(:,1),D(:,2),Energy,'nearest');
C1 = interp1(D(:,1),D(:,3),Energy,'nearest');
C2 = interp1(D(:,1),D(:,4),Energy,'nearest');

CrossSec = (C0 + C1.*Energy + C2.*Energy.^2).*Energy.^-3 .* 1e-24;   % [cm^2]

Tau  = CrossSec.*ColDen;
Tran = exp(-Tau);
