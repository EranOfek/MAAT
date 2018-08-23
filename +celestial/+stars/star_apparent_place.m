function [NewCat]=star_apparent_place(Catalog,Equinox,Epoch,ObsEpoch,GeodCoo)
%------------------------------------------------------------------------------
% star_apparent_place function                                           ephem
% Description: Compute the apparent place of stars at a given epoch,
%              given the star mean place, proper motion, parallax,
%              radial velocity, at a reference epoch T0.
% Input  : - Star catalog:
%            [RA, Dec, Mu_RA, Mu_Dec, Parallax, RV].
%            where RA, Dec are the equatorial coordinates of the stars in
%            radians.
%            Mu_RA, Mu_Dec are the proper motion in arcsec per julian century
%            (NaN if unknown).
%            Parallax is the star parallax in arcsec (NaN if unknown).
%            RV is the star radial velocity in km/sec (NaN if unknown).
%          - Catalog equinox (JD or 'J2000' | 'J1991.25' | 'B1950.0').
%          - Catalog epoch (JD or 'J2000' | 'J1991.25' | 'B1950.0')
%          - Epoch of observation: JD in the TT time scale or date
%            [D M Y frac_day]
%          - Observer geodetic coordinates: [Long, Lat, Height] in
%            (rad, rad, meters).
%            Default is [NaN, NaN, NaN] in this case diurnal aberration
%            is neglected.
% Output : - [R.A., Dec.] in radians referred to equinox and epoch of
%            observation.
% Reference: Explanatory supplement P. 121, 152
% Tested : Matlab 5.3
%     By : Eran O. Ofek                   October 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------------------
if (nargin==4),
   GeodCoo = [NaN, NaN, NaN];
elseif (nargin==5),
   % do nothing
else
   error('Illigal number of inout arguments');
end

RAD = 180./pi;
if (isstr(Equinox)==1),
   switch Equinox,
    case 'J2000.0'
       Equinox = 2451545.0;
    case 'J1991.25'
       Equinox = 
    case 'B1950.0'
       Equinox = 
    otherwise
       error('Unknown Equinox');
   end
end

if (isstr(Epoch)==1),
   switch Epoch,
    case 'J2000.0'
       Epoch = 2451545.0;
    case 'J1991.25'
       Epoch = 
    case 'B1950.0'
       Epoch = 
    otherwise
       error('Unknown Equinox');
   end
end

if (length(ObsEpoch)==4),
   ObsEpoch = julday(ObsEpoch);
end

SizeCat  = size(Catalog);
Nt       = SizeDate(1);

% prep catalog - convert to raw vectors
RA       = Catalog(:,1).';   % rad
Dec      = Catalog(:,2).';   % rad
Mu_RA    = Catalog(:,3).';   % "/cy
Mu_Dec   = Catalog(:,4).';   % "/cy
Parallax = Catalog(:,5).';   % "
RV       = Catalog(:,6).';   % km/sec

I            = find(isnan(Mu_RA)==1);
Mu_RA(I)     = 0;
I            = find(isnan(Mu_Dec)==1);
Mu_Dec(I)    = 0;
I            = find(isnan(Parallax)==1);
Parallax(I)  = 1e-7;  % 10 Mpc
I            = find(isnan(RV)==1);
RV(I)        = 0;



% convert JD_TT to JD_TDB if needed

% Find the rotation matrix : N(t)*P(t)
% for converting coordinates from mean equinox (Equinox) and
% equator (Equinox) to the epoch of observation (ObsEpoch).
if (Equinox==2451545.0),
   Rot = rotm_coo('Pd',ObsEpoch);
else
   Rot = rotm_coo('p',Equinox) * rotm_coo('Pd',ObsEpoch);
end

% Vector of stars mean  place at reference epoch (Epoch),
% represented as a 3-D vector in AU, with origin at the
% solar system Barycenter, and referred to the mean
% equator and equinox of catalog (Equinox).
% (One vector per raw).
Dist   = 1./sin(Parallax);  % (au) parallax should be in radians.
UB_T0  = [Dist.*cos(Dec).*cos(RA); Dist.*cos(Dec).*sin(RA); Dist.*sin(Dec)];

% find the space motion vector (in AU/day) of the star at the
% reference epoch (Epoch), obtained from the catalog proper motions,
% parallax, and radial velocity.
S      = 2.*pi./(360.*3600.*36525);  % convert from "/cy to radians/day
K      = 86400./(1.49597870660e8);   % convert from km/s to au/day
UtB_T0 = zeros(3,Nt);
UtB_T0(1,:) = -cos(Dec).*sin(RA).*15   .*S.*Dist    .*Mu_RA  ...
              -sin(Dec).*cos(RA)       .*S.*Dist    .*Mu_Dec ...
              +cos(Dec).*cos(RA)    .*K         .*RV;

UtB_T0(2,:) = +cos(Dec).*cos(RA).*15   .*S.*Dist    .*Mu_RA  ...
              -sin(Dec).*sin(RA)       .*S.*Dist    .*Mu_Dec ...
              +cos(Dec).*sin(RA)    .*K         .*RV;

UtB_T0(3,:) = +0                                             ...
              +cos(Dec)                .*S.*Dist    .*Mu_Dec ...
              +sin(Dec)             .*K         .*RV;
 

% find the Barycentric position of the Earth/Observer at the
% epoch of observation (ObsEpoch), referred to the
% mean equator and equinox (Equinox).
[EB, EBvel]    = calc_vsop87(ObsEpoch, 'Earth','e','E');  % J2000.0
% add topocentric velocity of observer (for diurnal aberration).
% Ignore DeltaT (should be: ObsEpoch - DeltaT)
if (isnan(GeodCoo(1))==0),
   [G,Gt]         = topocentric_vec(ObsEpoch,GeodCoo);
   PrecessV       = rotm_coo('pd',ObsEpoce)
   Gt             = PrecessV*Gt;
   EBvel          = EBvel + Gt;   % topocentric velocity
else
   % neglecting diurnal aberration
end

% precess to Equinox
if (Equinox==2451545.0),
   % already in J2000.0
else
   Precess = rotm_coo('P',Equinox);
   EB      = Precess*EB;
   EBvel   = Precess*EBvel;
end

P0 = UB_T0 + (JD_TDB - JD_Epoch).*UtB_T0;
P1 = UB_T0 + (JD_TDB - JD_Epoch).*UtB_T0 - EB;

% operate light deflection on P1:
ObsHelio = calc_vsop87(ObsEpoch,'Earth','a','E');   % Earth heliocentric position (J2000)
P2 = light_deflection(P1,P0,ObsHelio);

% operate aberration of light on P2:
P3 = light_abberation(P2,EBvel);

% convert to apparent place
P4 = zeros(3,Nt);
for I=1:1:Nt,
   P4(:,I) = Rot*P3(:,I);
end

% convert P4 to R.A. and Dec.
NewCat = cosined(P4.');


