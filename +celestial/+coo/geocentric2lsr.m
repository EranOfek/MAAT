function [Sp,Err,Gp]=geocentric2lsr(JD,Speed,SpeedErr,GalCoo,Type,VelLSR,ErrLSR)
% Geocentric or heliocentric velocity to velocity relative to the LSR
% Package: celestial.coo
% Description: Approximate conversion of geocentric or heliocentric
%              velocity to velocity relative to the local standard of
%              rest (LSR).
% Input  : - Column vector of JD.
%          - Column vector of measured speed
%            (+ for moving in observer direction).
%          - Column vector of error in measured speed.
%          - Two column matrix of Galactic coordinates [Long, Lat]
%            in radians.
%          - Conversion type:
%            'G' - Geocentric to LSR/Galactocentric.
%            'H' - Heliocentric to LSR/Galactocentric.
%          - The velocity of the Sun relative to the LSR
%            [u, v, w], in km/sec.
%            Default is [-9, 12, 7].
%          - Errors in the velocity of the Sun relative to the LSR.
%            Default is [0.5, 0.5, 0.5].
% Output : - Speed relative to the LSR.
%          - Error in speed relative to the LSR.
%          - Galactocentric speed.
% Reference: Allen's astrophysicl quantities
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Sep 2003
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [V,Ve,GV]=celestial.coo.geocentric2lsr(2451545,20,1,[1 1],'G')
% Reliable: 2
%------------------------------------------------------------------------------
DefVelLSR   = [-9,   12,   7];
DefErrLSR   = [0.5, 0.5, 0.5];
if (nargin==5)
    VelLSR = DefVelLSR;
    ErrLSR = DefErrLSR;
elseif (nargin==6)
    ErrLSR = DefErrLSR;
elseif (nargin==7)
    % do nothing
else
   error('Illegal number of input arguments');
end

[SunRA0,SunDec0,SunR0,SolarLong0] = celestial.SolarSys.suncoo(JD,'j');
[SunRA1,SunDec1,SunR1,SolarLong1] = celestial.SolarSys.suncoo(JD+0.001,'j');
SunGalCoo0                    = celestial.coo.coco([SunRA0,SunDec0],'j2000.0','g');
SunGalCoo1                    = celestial.coo.coco([SunRA1,SunDec1],'j2000.0','g');
switch Type
 case 'G'
    EarthSpeed                    = 29.5;   % [km/sec]
 case 'H'
    EarthSpeed                    = 0;   % [km/sec]
 otherwise
    error('Unknown Type Option');
end

SunDist0                      = celestial.coo.sphere_dist(SunGalCoo0(:,1),SunGalCoo0(:,2),GalCoo(:,1),GalCoo(:,2));
SunDist1                      = celestial.coo.sphere_dist(SunGalCoo1(:,1),SunGalCoo1(:,2),GalCoo(:,1),GalCoo(:,2));

if (SunDist1>SunDist0)
   Sign = -1;   % Earth moving away from object
else
   Sign = +1;
end

DeltaVgeo = Sign.*EarthSpeed.*sin(SunDist0);

DelU = VelLSR(:,1).*cos(-GalCoo(:,1)).*cos(GalCoo(:,2));
DelV = VelLSR(:,2).*sin(-GalCoo(:,1)).*cos(GalCoo(:,2));
DelW = VelLSR(:,3)                   .*sin(GalCoo(:,2));

Sp   = Speed + DelU + DelV + DelW + DeltaVgeo;
Err  = sqrt(SpeedErr.^2 + ErrLSR(:,1).^2 + ErrLSR(:,2).^2 + ErrLSR(:,3).^2);

Gp   = Speed + 9.*cos(GalCoo(:,1)).*cos(GalCoo(:,2)) + 232.*sin(GalCoo(:,1)).*cos(GalCoo(:,2)) + 7.*sin(GalCoo(:,2));

