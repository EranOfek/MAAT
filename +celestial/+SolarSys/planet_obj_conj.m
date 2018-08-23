function [Info]=planet_obj_conj(StartDate,EndDate,DistThresh,ObserverPlanet,Planet,ObjCoo)
%------------------------------------------------------------------------------------------
% planets_obj_conj function                                            ephem
% Description:   Calculate planet-object conjunctions/occultations,
%                                 In which two planets are occult or found within
%                                 DistThresh from each other as observed from
%                                 ObserverPlanet.
% Input  : - Start of search date or JD.
%            If one column vector is given then one JD.
%            If 4 columns matrix is given then [D M Y frac_day] per line.
%          - End of search date or JD.
%            If one column vector is given then one JD.
%            If 4 columns matrix is given then [D M Y frac_day] per line.
%          - Maximum angular distance threshold in arcsec, default is 3600.
%          - Observer's planet:
%            'Mercury' | 'Venus' | 'Earth' |  'Mars' |
%            'Jupiter' | 'Saturn' | 'Uranus' | 'Neptune'
%            Default is 'Earth'.
%          - Planet: 
%            'Mercury' | 'Venus' | 'Earth' |  'Mars' |
%            'Jupiter' | 'Saturn' | 'Uranus' | 'Neptune'
%          - Second object J2000.0 coordinates [RA, Dec] in radians.
% Output : - General info about conjunction/occultation
%            [JD_ExMin, Date, MinDist.*3600.*RAD, PA_MinDist.*RAD, RA, Dec, Elon.*RAD, Delta, HP.*RAD.*3600, Occ]
% Tested : Matlab 5.3
%     By : Eran O. Ofek          December 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------------------------------
AU        = 149597870.66; %km
RAD       = 180./pi;

DistThresh = DistThresh./3600./RAD;  % convert to radians

ObserverPlanet = 'Earth';

% safty angular distance: include, low accuracy error + planet size + parallax + time resolution
SaftyDist = 0.3./RAD;


InterpMethod = 'linear';

if (length(StartDate)==4),
   StartJD = julday([StartDate]);
elseif (length(StartDate)==1),
   StartJD = StartDate;
else
   error('StartDate should be in JD or Date format');
end

if (length(EndDate)==4),
   EndJD = julday([EndDate]);
elseif (length(EndDate)==1),
   EndJD = EndDate;
else
   error('EndDate should be in JD or Date format');
end


%---------------------------------------------------
% calculate low accuracy ephemeris for the planets
%---------------------------------------------------
% in interval (days).
Interval = 0.1;
JD_Range = [StartJD-2:Interval:EndJD+2].';
N        = length(JD_Range);

[Coo,Dist,Ang,Mag] = planet_ephem(JD_Range,Planet,ObserverPlanet,'J2000');
[AngDist,PA]       = sphere_dist(ObjCoo(1),ObjCoo(2),Coo(:,1),Coo(:,2));

[JD_Min, ExactMinY] = local_min([JD_Range, AngDist]);

K                   = find(ExactMinY<DistThresh);
Nmin                = length(K);


Info = zeros(Nmin,2);

for I=1:1:Nmin,
   Ind = K(I);

   JD_ExMin    = JD_Min(Ind);
   Date        = jd2date(JD_ExMin);
   Date        = [Date(1:3), convertdms(Date(4),'f','H')];

   MinDist     = ExactMinY(Ind);
   PA_MinDist  = interp4deg([JD_Range,PA],JD_ExMin);
   % PlanetI R.A.
   RA          = interp4deg([JD_Range, Coo(:,1)],JD_ExMin);
   % PlanetI Dec.
   Dec         = interp4deg([JD_Range, Coo(:,2)],JD_ExMin);
   % approximate elongation from the Sun
   Elon        = interp4deg([JD_Range, Ang(:,3)],JD_ExMin);
   % PlanetI-observer distance
   Delta       = interp4deg([JD_Range, Dist(:,1)],JD_ExMin);
   % PlanetI equatorial and polar angular semi-diameter
   PlLat       = 0;
   [Phys, Ang] = planet_radius(Planet,Delta,PlLat);
   % PlanetI horizontal parallax
   [ObserverRadius]=planet_radius(ObserverPlanet);
   HP          = asin(ObserverRadius(1)./(Delta.*AU));
   

   % check if occultation is visible from somewhere on ObserverPlanet
   Occ         = MinAngDist<(Ang(1) + HP);

   Info(I,1:2) = [JD_ExMin, Date, MinDist.*3600.*RAD, PA_MinDist.*RAD, RA, Dec, Elon.*RAD, Delta, HP.*RAD.*3600, Occ]; 

end

