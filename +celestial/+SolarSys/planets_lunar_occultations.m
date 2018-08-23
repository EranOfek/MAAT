function [OccList,PlanetNames]=planets_lunar_occultations(StartJD,EndJD,PlanetNames,GeodCoo,MinSep)
% Search and calculates for lunar occultations of the planets
% Package: celestial.SolarSys
% Description: Calculate local circumstences for lunar occultations of
%              Planets and asteroids.
%              Only events in which the planet is above the local horizon
%              will be selected.
% Input  : - Search start date [JD] or [D M Y].
%          - Search end date [JD] or [D M Y].
%          - Cell array of planets to search.
%            Alternatively a string indicating one of the following options:
%            'planets' - search all planets {'Mercury',..,'Neptune'}
%            'all'     - search all planets and 4 primary asteroids and 4
%                        Jovian satellites.
%          - Geodetic position [Long, Lat, Height],
%            where Long and Lat are in deg, and Height in km.
%          - Minimum seperation conjunctions to find [deg].
% Output : - Matrix of occultations with the following columns:
%            (1) - Disappereance (1) or Reapperaence (2).
%            (2) - Planet index in the planet names cell array (see next output argument).
%            (3) - JD for outer contact [UTC]
%            (4) - JD for inner contact [UTC]
%            (5) - Planet J2000 RA [deg]
%            (6) - Planet J2000 Dec [deg]
%            (7) - PA of occultations [deg] relative to Moon center
%            (8) - Planet Azimuth [deg]
%            (9) - Planet Altitude [deg]
%            (10)- Planet Magnitude
%            (11)- Planet angular diameter [arcsec]
%            (12)- Moon Illuminated fraction []
%            (13) - Solar Elongation [deg]
%          - Cell array of correspondinf planet names
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Dec 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: RAD = 180./pi;
%          [OccList,PlanetNames]=planets_lunar_occultations(2451545,2451545+100,{'Venus'},[35 32 0]./RAD,1);
% Reliable: 2
%-----------------------------------------------------------------------------

import celestial.SolarSys.*
import celestial.coo.*

InterpMethod = 'linear';



RAD = 180./pi;

if (isstr(PlanetNames)==1),
   switch lower(PlanetNames)
    case 'planets'
       PlanetNames = {'Mercury','Venus','Mars','Jupiter','Saturn','Uranus','Neptune'};
    case 'all'
       PlanetNames = {'Mercury','Venus','Mars','Jupiter','Saturn','Uranus','Neptune','Ceres','Pallas','Juno','Vesta','Io','Europa','Ganymede','Callisto'};
    otherwise
       error('Unknown Planet option');
   end
end


StepSize = 1;
StepUnit = 'h';
%--- get ephemerids for the Moon ---
Moon = get_horizons(StartJD,EndJD,'Moon','Geod',GeodCoo,'StepSize',StepSize,'StepUnit',StepUnit);

%--- get ephemerids for a planet ---
Npl = length(PlanetNames);

Found = 0;
OccList = zeros(0,13);
for Ipl=1:1:Npl,
   PlanetNames{Ipl}
   Planet = get_horizons(StartJD,EndJD,PlanetNames{Ipl},'Geod',GeodCoo,'StepSize',StepSize,'StepUnit',StepUnit);

   List.Planet = [Planet.J_RA, Planet.J_Dec, Planet.Mag, Planet.Illum, Planet.AngDiam, Planet.AngSOT.*Planet.AngSOTd];
   List.Moon   = [Moon.J_RA, Moon.J_Dec, Moon.Mag, Moon.Illum, Moon.AngDiam, Moon.AngSOT.*Moon.AngSOTd];
   AngFlag     = [1 1 0 0 0 1];
   Col.AngDiam = 5;
   [Conj,IndMinSep]=search_conj(Planet.JD,List.Planet,[],AngFlag,List.Moon,[],AngFlag,MinSep./RAD);

   %--- Look for all conjunctions within MinSep ---
   
   %--- Look for occultations ---
   if (isempty(Conj)==1),
      Nconj = 0;
   else
      Nconj = length(Conj.MinJD);
   end
   for Iconj=1:1:Nconj,
      %--- check if occultation ---
      if (Conj.MinJD~=0 & Conj.MinDist(Iconj).*RAD.*3600<=(0.5.*(Conj.List1(Iconj,Col.AngDiam)+Conj.List2(Iconj,Col.AngDiam)))),
         %--- occultation found ---
         Found = Found + 1;
         Fine.Moon = get_horizons(Conj.MinJD(Iconj)-1,Conj.MinJD(Iconj)+1,'Moon','Geod',GeodCoo,'StepSize',1,'StepUnit','m');
         Fine.Pl   = get_horizons(Conj.MinJD(Iconj)-1,Conj.MinJD(Iconj)+1,PlanetNames{Ipl},'Geod',GeodCoo,'StepSize',1,'StepUnit','m');

         Fine.Dist = sphere_dist(Fine.Moon.J_RA,Fine.Moon.J_Dec,Fine.Pl.J_RA,Fine.Pl.J_Dec);
         Fine.DistOuter = Fine.Dist-(Fine.Moon.AngDiam./2+Fine.Pl.AngDiam./2)./(RAD.*3600);
         Fine.DistInner = Fine.Dist-(Fine.Moon.AngDiam./2-Fine.Pl.AngDiam./2)./(RAD.*3600);
         ListZeros = find_local_zeros(Fine.Moon.JD,Fine.DistOuter);
         if (size(ListZeros,1)==2),
            TimesOuter = ListZeros(1:2,1);
         else
            TimesOuter = [NaN; NaN];
         end
         ListZeros = find_local_zeros(Fine.Moon.JD,Fine.DistInner);
         if (size(ListZeros,1)==2),
            TimesInner = ListZeros(1:2,1);
         else
            TimesInner = [NaN; NaN];
         end

         % positions at Outer contact:
         Pos.Moon.RA  =  interp1(Fine.Moon.JD,Fine.Moon.J_RA,TimesOuter,InterpMethod);
         Pos.Moon.Dec =  interp1(Fine.Moon.JD,Fine.Moon.J_Dec,TimesOuter,InterpMethod);
         Pos.Pl.RA    =  interp1(Fine.Pl.JD,Fine.Pl.J_RA,TimesOuter,InterpMethod);
         Pos.Pl.Dec   =  interp1(Fine.Pl.JD,Fine.Pl.J_Dec,TimesOuter,InterpMethod);
         Pos.Pl.Mag   =  interp1(Fine.Pl.JD,Fine.Pl.Mag,TimesOuter,InterpMethod);         
         Pos.Pl.AngDiam  =  interp1(Fine.Pl.JD,Fine.Pl.AngDiam,TimesOuter,InterpMethod);
         Pos.Moon.Illum  =  interp1(Fine.Moon.JD,Fine.Moon.Illum,TimesOuter,InterpMethod);

         [Pos.Dist,Pos.PA] = sphere_dist(Pos.Moon.RA,Pos.Moon.Dec,Pos.Pl.RA,Pos.Pl.Dec);

         [Pos.Sun.RA1, Pos.Sun.Dec1] = suncoo(TimesOuter(1),'j');
         [Pos.Sun.RA2, Pos.Sun.Dec2] = suncoo(TimesOuter(2),'j');
         [SunElong(1)] = sphere_dist(Pos.Sun.RA1, Pos.Sun.Dec1, Pos.Pl.RA(1), Pos.Pl.Dec(1)).*RAD;
         [SunElong(2)] = sphere_dist(Pos.Sun.RA2, Pos.Sun.Dec2, Pos.Pl.RA(2), Pos.Pl.Dec(2)).*RAD;

         HorizCoo = horiz_coo([Pos.Pl.RA, Pos.Pl.Dec],TimesOuter,GeodCoo(1:2)./RAD);
         % disapearance:
         L = 1;
         OccList(Found,:) = [1, Ipl, TimesOuter(L), TimesInner(L), ...
                             Pos.Pl.RA(L).*RAD, Pos.Pl.Dec(L).*RAD, Pos.PA(L).*RAD, HorizCoo(L,:).*RAD, ...
			      Pos.Pl.Mag(L), Pos.Pl.AngDiam(L), Pos.Moon.Illum(L), SunElong(L)];
         % Reapereance:
         Found = Found + 1;
         L = 2;
         OccList(Found,:) = [2, Ipl, TimesOuter(L), TimesInner(L), ...
                             Pos.Pl.RA(L).*RAD, Pos.Pl.Dec(L).*RAD, Pos.PA(L).*RAD, HorizCoo(L,:).*RAD, ...
			      Pos.Pl.Mag(L), Pos.Pl.AngDiam(L), Pos.Moon.Illum(L),SunElong(L)];

      end
   end
end
