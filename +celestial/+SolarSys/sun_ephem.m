function [Coo,Dist,Mag]=sun_ephem(JD_TT, ObserverType, Reference, ObsGeoCoo, DeltaT, Xp, Yp)
% Sun ephemeris
% POackage: celestial.SolarSys
% Description:  Sun ephemerids generator
% Input  : - Column vector of JD in the TT time scale, or date (see julday.m).
%          - Obsever position :
%            'Mercury' : Mercury
%            'Venus'   : Venus
%            'Earth'   : Earth, Geocentric (default)
%            'Topo'    : Earth, Topocentric   <-------- BUG
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
%          - Column vector of Delta T (=TT-UT1), needed for topocentric coordinates.
%          - The X coordinate of the celestial ephemeris pole
%            with respect to the terrestrial pole measured along
%            the 0 meridian (radians). - default is 0.
%            (one element per time).
%          - The Y coordinate of the celestial ephemeris pole
%            with respect to the terrestrial pole measured along
%            the 270 meridian (radians). default is 0.
%            (one element per time).
% Output : - Sun's coordinates:
%            [RA (rad), Dec (rad)]
%          - Sun's distances: [R (au), Tau (day)]
%            Tau is the observer-planet light time correction.
%          - [Mag], Sun apparent magnitude.
% Tested : Matlab 5.3 - in process
%     By : Eran O. Ofek                    May 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%----------------------------------------------------------------------------
SunAbsMag = -26.74;
AU        = 149597870.66; %km
DaySec    = 86400; % sec in day
C         = 173.1446327205364; % speed of light [au/day] 
MUC2      = 9.8706286e-9;  % mu/c^2 [au] 
RefEllips = 'WGS84';
LightTimeThreshold = 1e-9;

if (size(JD_TT,2)==1)
   JD_TT = JD_TT;
else
   JD_TT = celestial.timt.julday(JD_TT).';
end
N = length(JD_TT);
% convert TDB...



if (nargin==1)
   ObserverType = 'Earth';
   Reference    = 'J2000';
   ObsGeoCoo       = [NaN, NaN, NaN];
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==2)
   Reference    = 'J2000';
   ObsGeoCoo       = [NaN, NaN, NaN];
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==3)
   ObsGeoCoo       = [NaN, NaN, NaN];
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==4)
   DeltaT       = zeros(N,1);
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==5)
   Xp           = zeros(N,1);
   Yp           = zeros(N,1);
elseif (nargin==6)
   Yp           = zeros(N,1);
elseif (nargin==7)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (strcmp(ObserverType,'Topo')==1 && sum(isnan(ObsGeoCoo))>0)
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
    [ObsCoo,ObsVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Mercury', VSOP87_Type, VSOP87_OutCoo); 
 case 'Venus'
    [ObsCoo,ObsVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Venus', VSOP87_Type, VSOP87_OutCoo); 
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
        % do nothing - allready in date (Equatorial)
     case 'J2000'
        % rotate topocentric coordinates from Eq. date to Eq. J2000
        for I=1:1:N        
           Rot = celestial.coo.rotm_coo('pd',JD_TT(I));
           G(:,I)  = Rot*G(:,I);
           Gt(:,I) = Rot*Gt(:,I);
        end
     case 'ecliptic'
        error('Topocentric - ecliptic coordinates not available yet');
        % do nothing - allready in date (Ecliptic)
        % rotate topocentric coordinates from Eq. date to Ecl. date
        for I=1:1:N       
           Rot = celestial.coo.rotm_coo('e',JD_TT(I));
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
    [ObsCoo,ObsVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Mars', VSOP87_Type, VSOP87_OutCoo); 
 case 'Jupiter'
    [ObsCoo,ObsVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Jupiter', VSOP87_Type, VSOP87_OutCoo); 
 case 'Saturn'
    [ObsCoo,ObsVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Saturn', VSOP87_Type, VSOP87_OutCoo); 
 case 'Uranus'
    [ObsCoo,ObsVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Uranus', VSOP87_Type, VSOP87_OutCoo); 
 case 'Neptune'
    [ObsCoo,ObsVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Neptune', VSOP87_Type, VSOP87_OutCoo); 
 case 'Sun'
    ObsCoo = zeros(3,N);  
 otherwise
    error('Unknown object name');
end

% Sun Barycentric coordinates
[SunCoo,SunVel]=celestial.SolarSys.calc_vsop87(JD_TT, 'Sun', VSOP87_Type, VSOP87_OutCoo);


% solve for light time correction (Tau)
Tau     = zeros(N,1);  % light-time correction [sec]
TauLast = ones(N,1).*LightTimeThreshold.*100;

MaxDiff = max(abs(TauLast-Tau));
while (MaxDiff>LightTimeThreshold)
   % Sun Barycentric coordinates
   [SunCooTau,SunVelTau]=celestial.SolarSys.calc_vsop87(JD_TT-Tau, 'Sun', VSOP87_Type, VSOP87_OutCoo);

   P    = SunCooTau - ObsCoo;
   
   AbsP = sqrt(sum((P.^2),1));

   % Tau includes the effect of gravitational retardation due to the Sun
   TauLast = Tau; 
   Tau     = (AbsP + (2.*MUC2).*log((AbsP)./(AbsP)))./C;
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
    P1=P;
    % abberation of light
    P2=celestial.coo.light_abberation(P1,ObsVel);
    % Apply precession and nutation to the proper direction
    P3 = zeros(size(P2));
    for I=1:1:N
       R = celestial.coo.rotm_coo('Pd',JD_TT(I));
       P3(:,I) = R*P2(:,I);
    end
    % P3 is univt vector, apply distance
    P3 = P3.*[ones(3,1)*sqrt(sum(P.*P,1))];
 case 'ecliptic'
    % coordinates are allready in Ecliptic of date
    % light deflection
    P1=P;
    % abberation of light
    P2=celestial.coo.light_abberation(P1,ObsVel);
    P3 = P2.*[ones(3,1)*sqrt(sum(P.*P,1))];
    
 otherwise
    error('Unknown type of Reference system');
end


% --- Ephemerids ---
% planet equatorial coordinates
if (nargout>=1)
   Coo    = cosined(P3.');
end
if (nargout>=2)
   % Sun-observer distance
   S_RadVec = [sqrt(sum(P3.^2,1))].';
   
   Dist = [S_RadVec, Tau];

   Mag  = SunAbsMag + 5.*log10(S_RadVec);

end

