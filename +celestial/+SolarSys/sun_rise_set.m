function [Time,Ang]=sun_rise_set(Date,ObsCoo,TimeZone,DeltaT)
% Calculate Sun rise/set
% Package: celestial.SolarSys
% Description: Given the coordinates and observer position, calculate
%              rise/set/transit/twilight times and azimuth and altitude
%              for the Sun. The accuracy depends on the function used
%              for calculating the solar position. With the default
%              sun-position function, the geometric accuracy is about
%              a few seconds.
% Input  : - Date [D M Y] or julian day [JD], one date/jd per line.
%          - Observer geodetic position [East Long, North Lat, Height (meters)],
%            in radians.
%            Default is the Wise obs. position [  34.763, 30.596]/RAD.
%          - East Time Zone, default is 2 [hours].
%          - Delta T (=UT-UTC) [fraction of day], default is 0.
% Output : - [Morning Astronomical Twilight,
%             Morning Nautical Twilight,
%             Morning Civil Twilight,
%             Rise,
%             Transit,
%             Set,
%             Evening Civil Twilight,
%             Evening Nautical Twilight,
%             Evening Astronomical Twilight]
%            Line per each date. In fraction of day.
%            The time is relative to Local Time
%            (=UT + Time Zone) [hours].
%          - [Rise Az, Transit Alt, Set Az]
%            Line per each date. In radians.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Sep 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Time,Ang]=celestial.SolarSys.sun_rise_set(2451545+(0:5:365)',[35 32 0].*pi./180,0,0)
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;

if (nargin==1),
   ObsCoo = [  34.763, 30.596]./RAD;
   TimeZone = 2;
   DeltaT = 0;
elseif (nargin==2),
   TimeZone = 2;
   DeltaT = 0;
elseif (nargin==3),
   DeltaT = 0;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end

%DeltaT = 1./1440

SizeDate = size(Date);
if (SizeDate(2)==3),
   JD = celestial.time.julday([Date, zeros(SizeDate(1))]);
else
   JD = Date;
end

N = length(JD);

if (length(ObsCoo)>2),
   Height = ObsCoo(3);
else
   Height = 0;
end


RS_Alt  = (-0.8333 - 2.12./60.*sqrt(Height))./RAD;
T06_Alt = -6./RAD;
T12_Alt = -12./RAD;
T18_Alt = -18./RAD;

% Sun's apparent coordinates
[SunRA_p0, SunDec_p0] = celestial.SolarSys.suncoo(JD,'a');
[SunRA_m2, SunDec_m2] = celestial.SolarSys.suncoo(JD-2,'a');
[SunRA_m1, SunDec_m1] = celestial.SolarSys.suncoo(JD-1,'a');
[SunRA_p1, SunDec_p1] = celestial.SolarSys.suncoo(JD+1,'a');
[SunRA_p2, SunDec_p2] = celestial.SolarSys.suncoo(JD+2,'a');

I           = find(SunRA_p0<0);
SunRA_p0(I) = 2.*pi + SunRA_p0(I);
I           = find(SunRA_p1<0);
SunRA_p1(I) = 2.*pi + SunRA_p1(I);
I           = find(SunRA_p2<0);
SunRA_p2(I) = 2.*pi + SunRA_p2(I);
I           = find(SunRA_m1<0);
SunRA_m1(I) = 2.*pi + SunRA_m1(I);
I           = find(SunRA_m2<0);
SunRA_m2(I) = 2.*pi + SunRA_m2(I);

% Greenwich Apparent Sidereal Time at 0 TDT.
GAST  = celestial.time.lst(JD, 0, 'a');
GASTp = celestial.time.lst(JD+1, 0, 'a');
GASTm = celestial.time.lst(JD-1, 0, 'a');

Time = zeros(N,9);
Ang  = zeros(N,3);

for I=1:1:N,
   SunRA  = [SunRA_m1(I) ; SunRA_p0(I) ; SunRA_p1(I) ];
   if (SunRA(2)<SunRA(1)),
      SunRA(2) = SunRA(2) + 2.*pi;
   end
   if (SunRA(3)<SunRA(2)),
      SunRA(3) = SunRA(3) + 2.*pi;
   end
   SunDec = [SunDec_m1(I); SunDec_p0(I); SunDec_p1(I)];

   SunRA_m  = [SunRA_m2(I) ; SunRA_m1(I) ; SunRA_p0(I) ];
   if (SunRA_m(2)<SunRA_m(1)),
      SunRA_m(2) = SunRA_m(2) + 2.*pi;
   end
   if (SunRA_m(3)<SunRA_m(2)),
      SunRA_m(3) = SunRA_m(3) + 2.*pi;
   end
   SunDec_m = [SunDec_m2(I); SunDec_m1(I); SunDec_p0(I)];

   SunRA_p  = [SunRA_p0(I) ; SunRA_p1(I) ; SunRA_p2(I) ];
   if (SunRA_p(2)<SunRA_p(1)),
      SunRA_p(2) = SunRA_p(2) + 2.*pi;
   end
   if (SunRA_p(3)<SunRA_p(2)),
      SunRA_p(3) = SunRA_p(3) + 2.*pi;
   end
   SunDec_p = [SunDec_p0(I); SunDec_p1(I); SunDec_p2(I)];

   %--- Rise/Set
   [T_RS,Alt_RS,Az_RS] = celestial.SolarSys.rise_set([SunRA, SunDec], GAST(I) ,RS_Alt,  ObsCoo, DeltaT);
   T_RS = T_RS + TimeZone./24;
   if (T_RS(2)>1),
      [T_RSp,Alt_RSp,Az_RSp] = celestial.SolarSys.rise_set([SunRA_p, SunDec_p], GASTp(I) ,RS_Alt,  ObsCoo, DeltaT);
      T_RS(2)  = T_RSp(2);
      Az_RS(2) = Az_RSp(2);
   end
   if (T_RS(1)<0),
      [T_RSm,Alt_RSm,Az_RSm] = celestial.SolarSys.rise_set([SunRA_m, SunDec_m], GASTm(I) ,RS_Alt,  ObsCoo, DeltaT);
      T_RS(1) = T_RSm(1);
      Az_RS(1) = Az_RSm(1);
   end


   %--- Civil Twilight
   [T_06,Alt_06,Az_06] = celestial.SolarSys.rise_set([SunRA, SunDec], GAST(I) ,T06_Alt, ObsCoo, DeltaT);
   T_06 = T_06 + TimeZone./24;
   if (T_06(2)>1),
      [T_06p,Alt_06p,Az_06p] = celestial.SolarSys.rise_set([SunRA_p, SunDec_p], GASTp(I) ,T06_Alt,  ObsCoo, DeltaT);
      T_06(2) = T_06p(2);
   end
   if (T_06(1)<0),
      [T_06m,Alt_06m,Az_06m] = celestial.SolarSys.rise_set([SunRA_m, SunDec_m], GASTm(I) ,T06_Alt,  ObsCoo, DeltaT);
      T_06(1) = T_06m(1);
   end


   %--- Nautical Twilight
   [T_12,Alt_12,Az_12] = celestial.SolarSys.rise_set([SunRA, SunDec], GAST(I) ,T12_Alt, ObsCoo, DeltaT);
   T_12 = T_12 + TimeZone./24;
   if (T_12(2)>1),
      [T_12p,Alt_12p,Az_12p] = celestial.SolarSys.rise_set([SunRA_p, SunDec_p], GASTp(I) ,T12_Alt,  ObsCoo, DeltaT);
      T_12(2) = T_12p(2);
   end
   if (T_12(1)<0),
      [T_12m,Alt_12m,Az_12m] = celestial.SolarSys.rise_set([SunRA_m, SunDec_m], GASTm(I) ,T12_Alt,  ObsCoo, DeltaT);
      T_12(1) = T_12m(1);
   end



   %--- Astronomical Twilight
   [T_18,Alt_18,Az_18] = celestial.SolarSys.rise_set([SunRA, SunDec], GAST(I) ,T18_Alt, ObsCoo, DeltaT);
   T_18 = T_18 + TimeZone./24;
   if (T_18(2)>1),
      [T_18p,Alt_18p,Az_18p] = celestial.SolarSys.rise_set([SunRA_p, SunDec_p], GASTp(I) ,T18_Alt,  ObsCoo, DeltaT);
      T_18(2) = T_18p(2);
   end
   if (T_18(1)<0),
      [T_18m,Alt_18m,Az_18m] = celestial.SolarSys.rise_set([SunRA_m, SunDec_m], GASTm(I) ,T18_Alt,  ObsCoo, DeltaT);
      T_18(1) = T_18m(1);
   end

   Time(I,:) = [T_18(2), T_12(2), T_06(2), T_RS(2), T_RS(1), T_RS(3), T_06(3), T_12(3), T_18(3)];
   Ang(I,:)  = [Az_RS(1), Alt_RS, Az_RS(2)];
   

end


% debug printing
%convertdms(Time(1,:).','f','H')
