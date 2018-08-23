function [ConjList,PlanetNames]=plants_object_conjunctions(StartJD,EndJD,Object,PlanetNames,GeodCoo,MinSep)
% Local circumstences for conjunctions of planets with a given object
% Package: celestial.SolarSys
% Description: Calculate local circumstences for conjunctions of planets
%              with a given object.
% Input  : - Search start date [JD] or [D M Y].
%          - Search end date [JD] or [D M Y].
%          - J2000 [RA, Dec] of object in radians.
%          - Cell array of planets to search.
%            Alternatively a string indicating one of the following options:
%            'planets' - search all planets {'Mercury',..,'Neptune'}
%            'all'     - search all planets and 4 primary asteroids and 4
%                        Jovian satellites.
%          - Geodetic position [Long, Lat, Height],
%            where Long and Lat are in deg, and Height in km.
%          - Minimum seperation conjunctions to find [deg].
% Output : - Matrix of occultations with the following columns:
%            (1) - Planet index in the planet names cell array (see next output argument).
%            (2) - JD of conjunction [UTC]
%            (3) - Min distance [deg]
%            (4) - PA of conjunction [deg] relative to object
%            (5) - Solar-object angular distance [deg].
%          - Cell array of correspondinf planet names
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Dec 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/   
% Example: [ConjList,PlanetNames]=celestial.SolarSys.plants_object_conjunctions(2451545,2451545+10,[1 0],'planets',[1 1 0],1);
% Reliable: 2
%-----------------------------------------------------------------------------
InterpMethod = 'linear';


RAD = 180./pi;

if (ischar(PlanetNames)),
  switch lower(PlanetNames)
    case 'planets'
        PlanetNames = {'Mercury','Venus','Mars','Jupiter','Saturn','Uranus','Neptune'};
    case 'all'
        PlanetNames = {'Mercury','Venus','Mars','Jupiter','Saturn','Uranus','Neptune','Ceres','Pallas','Juno','Vesta','Io','Europa','Ganymede','Callisto'};
      otherwise
        error('Unknown Planet option');
   end
end
Npl = length(PlanetNames);

StepSize = 1;
StepUnit = 'h';

List.Obj    = [Object 0 0 0 0];
AngFlag     = [1 1 0 0 0 0];

Found = 0;
ConjList = zeros(0,5);
for Ipl=1:1:Npl,
   %PlanetNames{Ipl}
   Planet = celestial.SolarSys.get_horizons(StartJD,EndJD,PlanetNames{Ipl},'Geod',GeodCoo,'StepSize',StepSize,'StepUnit',StepUnit);

   List.Planet = [Planet.J_RA, Planet.J_Dec, Planet.Mag, Planet.Illum, Planet.AngDiam, Planet.AngSOT.*Planet.AngSOTd];

   [Conj,IndMinSep]=celestial.SolarSys.search_conj_sm(Planet.JD,List.Planet,[],AngFlag,List.Obj,[],AngFlag,MinSep./RAD);

   [Sun.RA, Sun.Dec] = celestial.SolarSys.suncoo(Conj.MinJD,'j');
   SunDist = celestial.coo.sphere_dist(Sun.RA,Sun.Dec, Object(1,1), Object(1,2));
   if (~isempty(Conj.MinJD)),
      N = length(Conj.MinJD);
      ConjList = [ConjList; [Ipl.*ones(N,1), Conj.MinJD, Conj.MinDist.*RAD, Conj.MinPA.*RAD, SunDist.*RAD]];
   end
end


I = (ConjList(:,3)<=MinSep & ConjList(:,2)>0);
ConjList = ConjList(I,:);
