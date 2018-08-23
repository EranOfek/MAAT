function [S,P,Table]=sfr(z)
% Estimate the cosmic star formation rate as a function of redshift
% Description: Fit the measured star formation rate as a function
%              of redshift and return its interpolated value in some
%              requested redshifts.
% Input  : - An array of redshifts at which to return the SFR.
% Output : - The star formation rate [SolarMass yr^-1 Mpc^-3]
%          - Best for polynomial coef.
%          - Table of measurments.
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [S,P,Table]=AstroUtil.cosmo.sfr((0:0.01:8)')
% Reliable: 2

FitOrder = 4;
% Visible light SFR estimation from:
% https://ned.ipac.caltech.edu/level5/March14/Madau/Madau5.html
Table.Cat = [0.01 0.1  -1.82 -0.02 +0.09
0.2  0.4  -1.50 -0.05 +0.05
0.4  0.6  -1.39 -0.08 +0.15
0.6  0.8  -1.20 -0.13 +0.31
0.8  1.2  -1.25 -0.13 +0.31
0.05 0.05 -1.77 -0.09 +0.08
0.05 0.2  -1.75 -0.18 +0.18
0.2  0.4  -1.55 -0.12 +0.12
0.4  0.6  -1.44 -0.10 +0.10
0.6  0.8  -1.24 -0.10 +0.10
0.8  1.0  -0.99 -0.08 +0.09
1.0  1.2  -0.94 -0.09 +0.09
1.2  1.7  -0.95 -0.08 +0.15
1.7  2.5  -0.75 -0.09 +0.49
2.5  3.5  -1.04 -0.15 +0.26
3.5  4.5  -1.69 -0.32 +0.22
0.92 1.33 -1.02 -0.12 +0.12
1.62 1.88 -0.75 -0.12 +0.12
2.08 2.37 -0.87 -0.09 +0.09
1.9  2.7  -0.75 -0.11 +0.09
2.7  3.4  -0.97 -0.15 +0.11
3.8  3.8  -1.29 -0.05 +0.05
4.9  4.9  -1.42 -0.06 +0.06
5.9  5.9  -1.65 -0.08 +0.08
7.0  7.0  -1.79 -0.10 +0.10
7.9  7.9  -2.09 -0.11 +0.11
7.0  7.0  -2.00 -0.11 +0.10
8.0  8.0  -2.21 -0.14 +0.14];
Col.zlow    = 1;
Col.zhigh   = 2;
Col.LogSFR  = 3;   % solarmass / yr / Mpc^3
Col.SFRerrL = 4;
Col.SFRerrH = 5;
Table.Col = Col;

P=polyfit((Table.Cat(:,1)+Table.Cat(:,2)).*0.5,Table.Cat(:,3),FitOrder);
S=10.^polyval(P,z);

