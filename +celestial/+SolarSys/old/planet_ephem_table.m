function planet_ephem_table(Year,ObjectName,ObserverType, Reference, ObsGeoCoo, TimeZone, OutTable,OutTableHTML);

% Input  : - Year 
%            Column vector of JD in the TT time scale.
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
%          - TimeZone in hours
%          - Output file name
%          - Output html file name

RAD    = 180./pi;
DefAlt = -0.5667/RAD;

JD_TT = julday([1 1 Year 0]) + [0:1:366].';
N = length(JD_TT);


% DeltaT in days
DeltaT = 68./86400;

FID     = fopen(OutTable,'w');
FIDhtml = fopen(OutTableHTML,'w');


fprintf(FIDhtml,'<html>');
fprintf(FIDhtml,'<head>');
fprintf(FIDhtml,'<title>%s ephemeids</title>',ObjectName);
fprintf(FIDhtml,'</head>');
fprintf(FIDhtml,'<body bgcolor=#ffffff text=#000000 link=#0000ff vlink=#ff0000>');

fprintf(FIDhtml,'<a href="http://astroclub.tau.ac.il/" target="_top">המועדון האסטרונומי</a> | <a href="http://astroclub.tau.ac.il" target="_top">אסטרופדיה</a>'); 
fprintf(FIDhtml,'<h3>%s ephemeids - %d </h3>',ObjectName,Year);
fprintf(FIDhtml,'<br>');
fprintf(FIDhtml,'<br>');
fprintf(FIDhtml,'<br>');
fprintf(FIDhtml,'<table>');
fprintf(FIDhtml,'<tr>')
fprintf(FIDhtml,'<td>תאריך');
fprintf(FIDhtml,'<td>עליה ישרה');
fprintf(FIDhtml,'<td>נטיה');
fprintf(FIDhtml,'<td>מרחק ארץ');
fprintf(FIDhtml,'<td>מרחק שמש');
fprintf(FIDhtml,'<td>אלונגציה');
fprintf(FIDhtml,'<td>בהירות');
fprintf(FIDhtml,'<td>פאזה');
fprintf(FIDhtml,'<td>קוטר זוויתי');
fprintf(FIDhtml,'<td>קבוצה');
fprintf(FIDhtml,'<td>זריחה');
fprintf(FIDhtml,'<td>שקיעה');
fprintf(FIDhtml,'<td>טרנזיט');


for I=1:1:N,
   %--- planet position ---
   [Coo,Dist,Ang,Mag] = planet_ephem(JD_TT(I)+0, ObjectName, ObserverType, Reference, ObsGeoCoo, DeltaT);
   [CooP]             = planet_ephem(JD_TT(I)-1, ObjectName, ObserverType, Reference, ObsGeoCoo, DeltaT);
   [CooN]             = planet_ephem(JD_TT(I)+1, ObjectName, ObserverType, Reference, ObsGeoCoo, DeltaT);

   %--- rise/set ---
   GAST               = lst(JD_TT(I),0);
   [TRS,TAlt,RS_Az]   = rise_set([CooP;Coo;CooN],GAST,DefAlt,ObsGeoCoo,DeltaT);
   [PosConst]         = constellation(Coo,JD_TT(I));
   Transit = convertdms(TRS(1)+TimeZone./24,'f','H');
   Rise    = convertdms(TRS(2)+TimeZone./24,'f','H');
   Set     = convertdms(TRS(3)+TimeZone./24,'f','H');

   Date = jd2date(JD_TT(I));

   %--- coordinated rounding ---
   Date = Date(1:3);
   RA   = convertdms(Coo(1),'r','H');
   Dec  = convertdms(Coo(2),'r','D');
   Dec(4) = round(Dec(4));
   if (Dec(1)==-1),
      DecS = '-';
   else
      DecS = '+';
   end

   %--- Angular diameter ---
   [PhysRad, AngRad] = planet_radius(ObjectName,Dist(1),0);

   %--- Print to file ---
   fprintf(FID,'%02d-%02d-%04d   %02d:%02d:%04.1f %c%02d:%02d:%02d  %6.3f %6.3f %6.1f  %5.1f  %5.2f  %2d  %c%c%c  %02d:%02d %02d:%02d %02d:%02d\n',Date, RA, DecS, Dec(2:4),Dist(1), Dist(2), Ang(3).*RAD, Mag(1), Mag(2), round(2.*AngRad(1).*RAD.*3600), PosConst{1}(1), PosConst{1}(2), PosConst{1}(3), Rise(1:2), Set(1:2), Transit(1:2) );

   fprintf(FIDhtml,'<tr>');
   fprintf(FIDhtml,'<td>%02d-%02d-%04d   <td>%02d:%02d:%04.1f <td>%c%02d:%02d:%02d  <td>%6.3f <td>%6.3f <td>%6.1f  <td>%5.1f  <td>%5.2f  <td>%2d  <td>%c%c%c  <td>%02d:%02d <td>%02d:%02d <td>%02d:%02d\n',Date, RA, DecS, Dec(2:4),Dist(1), Dist(2), Ang(3).*RAD, Mag(1), Mag(2), round(2.*AngRad(1).*RAD.*3600), PosConst{1}(1), PosConst{1}(2), PosConst{1}(3), Rise(1:2), Set(1:2), Transit(1:2) );
   fprintf(FIDhtml,'</tr>');

end

fprintf(FIDhtml,'</table>');
fprintf(FIDhtml,'<br>');
fprintf(FIDhtml,'</body>');
fprintf(FIDhtml,'</html>');

fclose(FID);
fclose(FIDhtml);
