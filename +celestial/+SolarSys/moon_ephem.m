function [Coo,Dist,Ang,K]=moon_ephem(JD_TT, Reference, ObserverType, ObsGeoCoo, DeltaT, Xp, Yp)
% ELP2000-82 Moon ephemeris
% Package: celestial.SolarSys
% Description: Calculate very accurate ELP2000-82 apparent coordinates of 
%              the moon as observed from one of the solar system planet's.
%              Assuming no light deflection.
% Input  : - Column vector of JD in the TT time scale.
%            (If needed in UTC time scale add DeltaT directly to JD).
%          - Reference frame:
%            'J2000'   : Equatorial, Mean Equinox and Equator of
%                        J2000.0 (FK5) (default).
%                        Not corrected for light deflection and abberation,
%                        but corrected for light-time and retardation of light.
%            'date'    : Equatorial, True Equinox and Equator of date
%            'ecliptic': Ecliptic, apparent True Equinox and ecliptic of date.
%          - Obsever position :
%            'Earth'   : Earth, Geocentric (default)
%            'Topo'    : Earth, Topocentric
%          - Observer WGS84 Geodetic position (for Topocentric coordinates).
%            [East Long. (rad), North Lat. (rad), Height (m)]
%          - Column vector of Delta T (=TT-UT1) in days,
%            needed for topocentric coordinates.
%          - The X coordinate of the celestial ephemeris pole
%            with respect to the terrestrial pole measured along
%            the 0 meridian (radians). - default is 0.
%            (one element per time).
%          - The Y coordinate of the celestial ephemeris pole
%            with respect to the terrestrial pole measured along
%            the 270 meridian (radians). default is 0.
%            (one element per time).
% Output : - Moon's coordinates:
%            [RA (rad), Dec (rad)]
%          - Moon's distances:
%            [Delta (au), R (au), obs-sun (au), Tau (day)]
%            where R is the Sun-Planet distance
%            and Delta is the Planet-Observer distance.
%            Tau is the observer-planet light time correction.
%          - Moon's angles (approximate):
%            [Theta (rad), Phi (rad), Elon (rad)]
%            Theta is the Observer-Sun-Moon angle.
%            Phi is the observer-planet-sun angle.
%            Elon is the sun-observer-planet angle.
%          - [K], Moon's illuminated fraction of the disk.
% Reference : ELP2000-82B (Chapront-Touze & Chapront 1982).
% Tested : Matlab 5.3
%     By : Eran O. Ofek                   Sep 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: UNKNOWN
%---------------------------------------------------------------------------

AU        = 149597870.66; %km
DaySec    = 86400; % sec in day
C         = 173.1446327205364; % speed of light [au/day] 
MUC2      = 9.8706286e-9;  % mu/c^2 [au] 
RefEllips = 'WGS84';
LightTimeThreshold = 1e-9;

N = length(JD_TT);

if (nargin==1),
   Reference    = 'J2000';
   ObserverType = 'Earth';
   ObsGeoCoo    = [NaN, NaN, NaN];
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==2),
   ObserverType = 'Earth';
   ObsGeoCoo    = [NaN, NaN, NaN];
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==3),
   ObsGeoCoo    = [NaN, NaN, NaN];
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==4),
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==5),
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==6),
   Yp           = zeros(N,1);
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (strcmp(ObserverType,'Topo')==1 && sum(isnan(ObsGeoCoo))>0),
   error('Unknown observer topocentric coordinates');
end

VSOP87_Type   = 'e';
switch Reference
 case 'J2000'
    VSOP87_OutCoo = 'E';
 case 'date'
    VSOP87_OutCoo = 'E';
 case 'ecliptic'
    VSOP87_OutCoo = 'd'; % default coordinates    
 otherwise
    error('Unknown Reference type');
end


%--- Observer Position ---
switch ObserverType
 case 'Earth'
    [ObsCoo,ObsVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Earth', VSOP87_Type, VSOP87_OutCoo); 
 case 'Topo'
    [ObsCoo,ObsVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Earth', VSOP87_Type, VSOP87_OutCoo); 
    % convert geocentric to topocentric
    %[Geoc,GeocCart]=geod2geoc(ObsGeoCoo,RefEllips);
    JD_UT1 = JD_TT - DeltaT;
    [G,Gt]=celestial.coo.topocentric_vec(JD_UT1, ObsGeoCoo, RefEllips, Xp, Yp);

    switch Reference
     case 'date'
        % G,Gt - allready in true equator and equinox of date (Equatorial)
     case 'J2000'
        % rotate topocentric coordinates from Eq. date to Eq. J2000
        for I=1:1:N,        
           Rot = celestial.coo.rotm_coo('pd',JD_TT(I));
           G(:,I)  = Rot*G(:,I);
           Gt(:,I) = Rot*Gt(:,I);
        end
     case 'ecliptic'
        error('Topocentric - ecliptic coordinates not available yet');
        % do nothing - allready in date (Ecliptic)
        % rotate topocentric coordinates from Eq. date to Ecl. date
        for I=1:1:N,        
           Rot = celestial.coo.rotm_coo('e',JD_TT(I));
           G(:,I)  = Rot*G(:,I);
           Gt(:,I) = Rot*Gt(:,I);
        end
     otherwise
        error('Unknown Reference coordinates type');
    end

    ObsCoo = ObsCoo + G./(AU.*1000);
    ObsVel =          Gt.*DaySec/(AU.*1000);
 otherwise
    error('Unknown observer type');
end

Tau     = zeros(N,1);  % light-time correction [sec]
TauLast = ones(N,1).*LightTimeThreshold.*100;

MaxDiff = max(abs(TauLast-Tau));
while (MaxDiff>LightTimeThreshold),
   %--- Moon Position ---
   ELP82_Type  = 'q2000';
   [X,Y,Z]      = celestial.SolarSys.moonpos(JD_TT-Tau,ELP82_Type);
   ObjCooG = [X;Y;Z]./AU;  % convert km to au (relative to Earth)

   switch ObserverType
    case 'Earth'
       % Geocentric coordinates
       ObjCooT = ObjCooG;
    case 'Topo'
       % Topocentric coordinates
       ObjCooT = ObjCooG - G./(AU.*1000);
    otherwise
       error('Unknown observer type');
   end

   MoonObsDist = sqrt(sum(ObjCooT.^2,1));  % au
   TauLast     = Tau;
   Tau         = [MoonObsDist./C].';  % light time correction
   MaxDiff = max(abs(TauLast-Tau));
end



%--- Sun Barycentric coordinates ---
[SunCoo,SunVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Sun', VSOP87_Type, VSOP87_OutCoo);
% Sun Barycentric coordinates
[SunCooTau,SunVelTau]=celestial.SolarSys.calc_vsop87(JD_TT-Tau, 'Sun', VSOP87_Type, VSOP87_OutCoo);
   
P    = ObjCooT;    % moon in respect to Earth
E    = ObsCoo - SunCoo;   % Earth respect to the sun
Q    = ObjCooT + ObsCoo - SunCooTau;  % Moon respect to the sun

AbsP = sqrt(sum((P.^2),1));
AbsE = sqrt(sum((E.^2),1));
AbsQ = sqrt(sum((Q.^2),1));


switch Reference
 case 'J2000'
    % coordinates are already in Equatorial J2000.0 FK5
    %P1 = P;  % no light deflection
    %P2 = light_abberation(P1,ObsVel);
    P3 = P;

 case 'date'
    % convert to Equinox of date
    % no light deflection
    P1 = P;
    % abberation of light
    ObsVel = Gt;
    P2     = P1; %light_abberation(P1,ObsVel);
    % Apply precession and nutation to the proper direction
    P3 = zeros(size(P2));
    for I=1:1:N,
       R = rotm_coo('Pd',JD_TT(I));
       P3(:,I) = R*P2(:,I);
    end
    % P3 is univt vector, apply distance
    P3 = P3.*[ones(3,1)*sqrt(sum(P.*P,1))];
 case 'ecliptic'
    % coordinates are allready in Ecliptic of date
    % no light deflection
    P1 = P;
    % abberation of light
    ObsVel = Gt;
    P2 = P1; %light_abberation(P1,ObsVel);
    P3 = P2.*[ones(3,1)*sqrt(sum(P.*P,1))];
    
 otherwise
    error('Unknown type of Reference system');
end

% --- Ephemerids ---
% planet equatorial coordinates
if (nargout>=1),
   Coo    = cosined(P3.');
end

%--- Distances ---
if (nargout>=2)
   % planet-observer distance
   P_RadVec = [sqrt(sum(P3.^2,1))].';
   % planet-sun distance
   R_RadVec = [sqrt(sum(Q.^2,1))].';
   % observer radius vector
   O_RadVec = [sqrt(sum(E.^2,1))].';
   
   Dist = [P_RadVec, R_RadVec, O_RadVec, Tau];
end

%--- Angles ---
if (nargout>=3),
   % observer-planet-sun angle
   Phi   = acos(sum(P.*Q,1)./(AbsP.*AbsQ));

   % sun-observer-planet angle
   Elon = asin(AbsQ.*sin(Phi)./AbsE);

   % observer-sun-planet angle
   Theta = pi - Phi - Elon;
   
   Phi   = Phi.'; 
   Elon  = Elon.';
   Theta = Theta.';

   Theta = pi - Phi - Elon;
   Ang   = [Theta, Phi, Elon];
end

%--- Illuminated Fraction ---
if (nargout>=4),
   % Illuminated fraction of the disk of a planet
   K = 0.5.*(1 + cos(Phi));
end



