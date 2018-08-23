function [DeltaV,D,ObjMoonDist,K]=moon_sky_brightness(Date,ObjCoo,GeodPos,C_Ext,Vsky)
% Krisciunas & Schaefer (1991) sky brightness model due to the Moon
% Package: celestial.SolarSys
% Description: Given the date, object equatorial coordinates, and 
%              observer geodetic position, calculate the excess in
%              sky brightness (V-band) in the object celestial position.
%              The function utilize the algorithm by 
%              Krisciunas & Schaefer (1991).
% Input  : - Date [day, month, year, frac_day] or JD.
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
% See also: moon_sky_brightness1.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DeltaV,D,ObjMoonDist,K]=celestial.SolarSys.moon_sky_brightness(2451545,[1 1],[1 1],0.2,21)
% Reliable: 2
%------------------------------------------------------------------------------
RAD = 180./pi;



which geod2geoc
if (nargin==2),
   GeodPos = [34.7630./RAD, 30.5960./RAD, 0];
   C_Ext      = 0.3;
   Vsky       = 21.7;
elseif (nargin==3),
   C_Ext      = 0.3;
   Vsky       = 21.7;
elseif (nargin==4),
   Vsky       = 21.7;
elseif (nargin==5),
   % do nothing
else
   error('Illigal number of input arguments');
end


% Julian Day
if (length(Date(1,:))>1),
   JD = celestial.time.julday([Date(1),Date(2),Date(3),Date(4)]);
else
   JD = Date;
end

% Geodetic to Geocentric
[GeocPos]=celestial.Earth.geod2geoc(GeodPos,'WGS84');

% Moon & Sun Position
[MoonRA,MoonDec,MoonHP] = celestial.SolarSys.mooncool(JD,GeocPos);
[SunRA,SunDec,SunR]     = celestial.SolarSys.suncoo(JD,'a');
MoonR = asin(MoonHP).*6378.137./149597870.0; %AU

% Moon Elongation
D = celestial.coo.sphere_dist(MoonDec,SunDec,MoonRA,SunRA);
Alpha = RAD.*(pi - D);


% Selenographic elongation of the Earth from the Sun
I = atan2((SunR.*sin(D)),(MoonR - SunR.*cos(D)));

% Illuminated fraction
K = 0.5.*(1 + cos(I));




% convert obj coo. to horiz coo.
ObjHor = celestial.coo.horiz_coo(ObjCoo,JD,GeodPos,'h');
%ObjHor = ObjCoo;
Z = pi./2 - ObjHor(:,2);

% Moon horizontal coo.
MoonHor = celestial.coo.horiz_coo([MoonRA,MoonDec],JD,GeodPos,'h');
Z_Moon = pi./2 - MoonHor(:,2);


% Object-Moon distance
%ObjMoonDist = sphere_dist(ObjCoo(2),MoonDec,ObjCoo(1),MoonRA);
ObjMoonDist = celestial.coo.sphere_dist(ObjHor(:,2),MoonHor(:,2),ObjHor(:,1),MoonHor(:,1));



% Moon Illumination
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








