function [DeltaV,D,ObjMoonDist,K]=moon_sky_brightness_h(D,MoonHor,ObjHor)
%------------------------------------------------------------------------------
% moon_sky_brightness_h function                                         ephem
% Description: Given the horizontal coordinates of the Moon and an object
%              and the observer geodetic position, calculate the excess
%              in sky brightness (V-band) in the object celestial 
%              position. The function utilize the algorithm by
%              Krisciunas & Schaefer (1991).
% Input  : - Moon elongation [rad].
%          - Moon horizontal coordinates [Az, Alt] in radians.
%          - Object horizontal coordinates [Az, Alt] in radians.
%            if one element is given than assumed to be JD.
%          - Object apparent equatorial
%            coordinates [RA, Dec] in radians.
%          - Observer geodetic position [East_Long, Lat, Height],
%            radians and meters above ref ellips.
%            default is wise observatory position.
%          - Extinction coef. in V. (default is 0.3mag/airmass).
%          - Sky brightness in V. (default is 21.7 mag/sq. arcsec.).
% Output : - The change in the V-band sky brightness caused by moonlight.
%          - Moon elongation, radians.
%          - Object-Moon distance, radians.
%          - Moon illuminated fraction.
% Reference : Krisciunas, K. and Schaefer, B. 1991 PASP 103, 1033.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    August 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: moon_sky_brightness.m
% Reliable: 2
%------------------------------------------------------------------------------
RAD = 180./pi;

%D        = [1];     % moon elongation
MoonDist = 384400;
%ObjHor   = [1 1];
%MoonHor  = [0 1];
C_Ext      = 0.3;
Vsky       = 21.7;



SunR = 1;

MoonR = MoonDist./150e6;

% Moon Elongation
Alpha = RAD.*(pi - D);


% Selenographic elongation of the Earth from the Sun
I = atan2((SunR.*sin(D)),(MoonR - SunR.*cos(D)));

% Illuminated fraction
K = 0.5.*(1 + cos(I));


Z = pi./2 - ObjHor(2);

% Moon horizontal coo.
Z_Moon = pi./2 - MoonHor(:,2);


% Object-Moon distance
%ObjMoonDist = distsp(ObjCoo(2),MoonDec,ObjCoo(1),MoonRA);
ObjMoonDist = distsp(ObjHor(2),MoonHor(2),ObjHor(1),MoonHor(1));



% moon ilumination
I_star = 10.^(-0.4.*(3.84 + 0.026.*abs(Alpha) + 4e-9.*Alpha.^4));
F_Rho  = (10.^5.36).*(1.06 + cos(ObjMoonDist).^2) + 10.^(6.15 - RAD.*ObjMoonDist./40);

Xz     = (1 - 0.96.*sin(Z).^2).^(-0.5);
XzMoon = (1 - 0.96.*sin(Z).^2).^(-0.5);

% moon sky brigtness in nanoLamberts
%Bmoon  = F_Rho.*I_star.*(1 - 10.^(-0.4.*C_Ext.*Xz)).*10.^(-0.4.*C_Ext.*Xz);
Bmoon  = F_Rho.*I_star.*(1 - 10.^(-0.4.*C_Ext.*Xz)).*10.^(-0.4.*C_Ext.*XzMoon);
IbelowHor = find(MoonHor(:,2)<0);
Bmoon(IbelowHor) = 0;

% convert nanLamberts to mag/sq. arcsec.
%Bsky   = -(log(Bmoon./34.08) - 20.7233)./0.92104


% convert sky brightness to nanoLamberts
VskyNL = 34.08.*exp(20.7233 - 0.92104.*Vsky);

DeltaV = -2.5.*log10((Bmoon + VskyNL)./VskyNL);









% distance function
function Dist=distsp(D1,D2,R1,R2)
Dist = acos(sin(D1).*sin(D2) + cos(D1).*cos(D2).*cos(R1-R2));
