function [Coo,Dist,Ang,Mag]=planet_ephem(JD_TT, ObjectName, ObserverType, Reference, ObsGeoCoo, DeltaT, Xp, Yp)
%----------------------------------------------------------------------------
% planet_ephem function                                                ephem
% Description: Planetary ephemerids generator based on VSOP87.
% Input  : - Column vector of JD in the TT time scale, or date.
%            (If needed in UTC time scale add DeltaT directly to JD).
%          - Object name:
%            'Mercury' : Mercury
%            'Venus'   : Venus
%            'Earth'   : Earth
%            'Mars'    : Mars
%            'Jupiter' : Jupiter
%            'Saturn'  : Saturn
%            'Uranus'  : Uranus
%            'Neptune' : Neptune
%          - Obsever position :
%            'Mercury' : Mercury
%            'Venus'   : Venus
%            'Earth'   : Earth, Geocentric (default)
%            'Topo'    : Earth, Topocentric
%            'Mars'    : Mars
%            'Jupiter' : Jupiter
%            'Saturn'  : Saturn
%            'Uranus'  : Uranus
%            'Neptune' : Neptune
%            'Sun'     : Sun (Heliocentric)
%          - Reference frame:
%            'J2000'   : Equatorial, Mean Equinox and Equator of
%                        J2000.0 (FK5) (default).
%                        Not corrected for light deflection and abberation,
%                        but corrected for light-time and retardation of light.
%            'date'    : Equatorial, True Equinox and Equator of date
%            'ecliptic': Ecliptic, apparent True Equinox and ecliptic of date.
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
% Output : - Planet's coordinates:
%            [RA (rad), Dec (rad)]
%          - Planet's distances:
%            [Delta (au), R (au), obs-sun (au), Tau (day)]
%            where R is the Sun-Planet distance
%            and Delta is the Planet-Observer distance.
%            Tau is the observer-planet light time correction.
%          - Planet's angles (approximate):
%            [Theta (rad), Phi (rad), Elon (rad)]
%            Theta is the Observer-Sun-Planet angle.
%            Phi is the observer-planet-sun angle.
%            Elon is the sun-observer-planet angle
%            The elongation is measured eastward.
%          - [Mag, K], Planet's magnitude and
%            Illuminated fraction of the disk of a planet.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    May 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Coo,Dist,Ang,Mag]=planet_ephem(2451545+[0:1:10].',...
%                                          'Venus',...
%                                          'Earth',...
%                                          'J2000');
%----------------------------------------------------------------------------
AU        = 149597870.66; %km
DaySec    = 86400; % sec in day
C         = 173.1446327205364; % speed of light [au/day] 
MUC2      = 9.8706286e-9;  % mu/c^2 [au] 
RefEllips = 'WGS84';
LightTimeThreshold = 1e-9;

if (size(JD_TT,2)==1),
   % do nothing
else
   JD_TT = julday(JD_TT).';
end
N = length(JD_TT);
% convert TDB...

if (nargin==2),
   ObserverType = 'Earth';
   Reference    = 'J2000';
   ObsGeoCoo    = [NaN, NaN, NaN];
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==3),
   Reference    = 'J2000';
   ObsGeoCoo    = [NaN, NaN, NaN];
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==4),
   ObsGeoCoo    = [NaN, NaN, NaN];
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==5),
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==6),
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==7),
   Yp           = zeros(N,1);
elseif (nargin==8),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (strcmp(ObserverType,'Topo')==1 & sum(isnan(ObsGeoCoo))>0),
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

ELP82_Type    = 'e2000';


switch ObserverType
 case 'Mercury'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Mercury', VSOP87_Type, VSOP87_OutCoo); 
 case 'Venus'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Venus', VSOP87_Type, VSOP87_OutCoo); 
 case 'Earth'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Earth', VSOP87_Type, VSOP87_OutCoo); 
 case 'Topo'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Earth', VSOP87_Type, VSOP87_OutCoo); 
    % convert geocentric to topocentric
    %[Geoc,GeocCart]=geod2geoc(ObsGeoCoo,RefEllips);
    JD_UT1 = JD_TT - DeltaT;
    [G,Gt]=topocentric_vec(JD_UT1, ObsGeoCoo, RefEllips, Xp, Yp);
    switch Reference
     case 'date'
        % do nothing - allready in date (Equatorial)
     case 'J2000'
        % rotate topocentric coordinates from Eq. date to Eq. J2000
        for I=1:1:N,        
           Rot = rotm_coo('pd',JD_TT(I));
           G(:,I)  = Rot*G(:,I);
           Gt(:,I) = Rot*Gt(:,I);
        end
     case 'ecliptic'
        error('Topocentric - ecliptic coordinates not available yet');
        % do nothing - allready in date (Ecliptic)
        % rotate topocentric coordinates from Eq. date to Ecl. date
        for I=1:1:N,        
           Rot = rotm_coo('e',JD_TT(I));
           G(:,I)  = Rot*G(:,I);
           Gt(:,I) = Rot*Gt(:,I);
        end        
     otherwise
        error('Unknown Reference coordinates type');
    end
    ObsCoo = ObsCoo + G./(AU.*1000);
    ObsVel = ObsVel + Gt.*DaySec/(AU.*1000);
 %case 'Moon'
 %   [X,Y,Z]=moonpos(JD_TT, ELP82_Type);
 %   ObsCoo = [X.';Y.';Z.'];
 case 'Mars'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Mars', VSOP87_Type, VSOP87_OutCoo); 
 case 'Jupiter'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Jupiter', VSOP87_Type, VSOP87_OutCoo); 
 case 'Saturn'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Saturn', VSOP87_Type, VSOP87_OutCoo); 
 case 'Uranus'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Uranus', VSOP87_Type, VSOP87_OutCoo); 
 case 'Neptune'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Neptune', VSOP87_Type, VSOP87_OutCoo); 
 case 'Sun'
    [ObsCoo,ObsVel]=calc_vsop87(JD_TT, 'Sun', VSOP87_Type, VSOP87_OutCoo);
    %ObsCoo = zeros(3,N);  
 otherwise
    error('Unknown object name');
end


% Sun Barycentric coordinates
[SunCoo,SunVel]=calc_vsop87(JD_TT, 'Sun', VSOP87_Type, VSOP87_OutCoo);

% solve for light time correction (Tau)
Tau     = zeros(N,1);  % light-time correction [sec]
TauLast = ones(N,1).*LightTimeThreshold.*100;

MaxDiff = max(abs(TauLast-Tau));
while (MaxDiff>LightTimeThreshold)
   switch ObjectName
    case 'Mercury'
       [ObjCoo,ObjVel]=calc_vsop87(JD_TT-Tau, 'Mercury', VSOP87_Type, VSOP87_OutCoo); 
    case 'Venus'
       [ObjCoo,ObjVel]=calc_vsop87(JD_TT-Tau, 'Venus', VSOP87_Type, VSOP87_OutCoo); 
    case 'Earth'
       [ObjCoo,ObjVel]=calc_vsop87(JD_TT-Tau, 'Earth', VSOP87_Type, VSOP87_OutCoo); 
    %case 'Moon'
    %   [X,Y,Z]=moonpos(JD_TT, ELP82_Type);
    %   ObjCoo = [X.';Y.';Z.'];
    case 'Mars'
       [ObjCoo,ObjVel]=calc_vsop87(JD_TT-Tau, 'Mars', VSOP87_Type, VSOP87_OutCoo); 
    case 'Jupiter'
       [ObjCoo,ObjVel]=calc_vsop87(JD_TT-Tau, 'Jupiter', VSOP87_Type, VSOP87_OutCoo); 
    case 'Saturn'
       [ObjCoo,ObjVel]=calc_vsop87(JD_TT-Tau, 'Saturn', VSOP87_Type, VSOP87_OutCoo); 
    case 'Uranus'
       [ObjCoo,ObjVel]=calc_vsop87(JD_TT-Tau, 'Uranus', VSOP87_Type, VSOP87_OutCoo); 
    case 'Neptune'
       [ObjCoo,ObjVel]=calc_vsop87(JD_TT-Tau, 'Neptune', VSOP87_Type, VSOP87_OutCoo);      
    otherwise
       error('Unknown object name');
   end
   
   
   % Sun Barycentric coordinates
   [SunCooTau,SunVelTau]=calc_vsop87(JD_TT-Tau, 'Sun', VSOP87_Type, VSOP87_OutCoo);
   
   P    = ObjCoo - ObsCoo;
   E    = ObsCoo - SunCoo;
   Q    = ObjCoo - SunCooTau;

   AbsP = sqrt(sum((P.^2),1));
   AbsE = sqrt(sum((E.^2),1));
   AbsQ = sqrt(sum((Q.^2),1));

   % Tau includes the effect of gravitational retardation due to the Sun
   TauLast = Tau; 
   Tau     = (AbsP + (2.*MUC2).*log((AbsE + AbsP + AbsQ)./(AbsE - AbsP + AbsQ)))./C;
   Tau     = Tau.';
   MaxDiff = max(abs(TauLast-Tau));
end

switch Reference
 case 'J2000'
    % coordinates are already in Equatorial J2000.0 FK5
    P3 = P;

 case 'date'
    % convert to Equinox of date
    % light deflection
    P1 = light_deflection(P,Q,E);
    % abberation of light
    P2 = light_abberation(P1,ObsVel);
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
    % light deflection
    P1 = light_deflection(P,Q,E);
    % abberation of light
    P2 = light_abberation(P1,ObsVel);
    P3 = P2.*[ones(3,1)*sqrt(sum(P.*P,1))];
    
 otherwise
    error('Unknown type of Reference system');
end

% --- Ephemerids ---
% planet equatorial coordinates
if (nargout>=1),
   Coo    = cosined(P3.');
end
if (nargout>=2)
   % planet-observer distance
   P_RadVec = [sqrt(sum(P3.^2,1))].';
   % planet-sun distance
   R_RadVec = [sqrt(sum(Q.^2,1))].';
   % observer radius vector
   O_RadVec = [sqrt(sum(E.^2,1))].';
   
   Dist = [P_RadVec, R_RadVec, O_RadVec, Tau];
end
if (nargout>=3),
   % observer-planet-sun angle
   Phi   = acos(sum(P.*Q,1)./(AbsP.*AbsQ));

   % sun-observer-planet angle
   %Elon = asin(AbsQ.*sin(Phi)./AbsE);
   Elon = acos((AbsP.^2 + AbsE.^2 - AbsQ.^2)./(2.*AbsP.*AbsE));
   % sign of elongation
   PlCooEq  = cosined(P.');
   SunCooEq = cosined(-E.');
   [DistPS,PA_PS] = sphere_dist(PlCooEq(:,1),PlCooEq(:,2),SunCooEq(:,1),SunCooEq(:,2));
   SignElon = sign(PA_PS - pi).';

   Elon = Elon.*SignElon;

   % observer-sun-planet angle
   Theta  = pi - Phi - abs(Elon);
   
   Phi   = Phi.';
   Elon  = Elon.';
   Theta = Theta.';
   Ang   = [Theta, Phi, Elon];
end
if (nargout>=4),
   Band = 'V';
   switch ObjectName
    case 'Saturn'
       [SLO,SLS,DU,P]=saturn_rings(JD_TT);
       % bug
       DU  = 0;
       SLO = 0;       
    otherwise
       DU  = 0;
       SLO = 0;
   end
   Mag = planets_magnitude(ObjectName, R_RadVec, P_RadVec, O_RadVec, Band, DU, SLO);
   % Illuminated fraction of the disk of a planet
   K = 0.5.*(1 + cos(Phi));
   Mag = [Mag, K];
end


