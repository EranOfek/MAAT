function Res=daily_observability(ObsCoo, Date, RA, Dec, Equinox)
% Observability plot fot a specific coordinates
% Package: telescope.obs
% Description: Plot the observability of a given object from a give location
%              on Earth during one night.
%              This program will plot the object Alt during the night,
%              along with the Sun/Moon alt and the excess in sky
%              brightness due to the Moon.
% Input  : - Observer Geodetic position, [E-Long, N-Lat] in radians.
%            or give string to choose from list (see observatory_coo.m):
%            'Wise' - wise observatory        [  34.763, 30.596]/RAD.
%            'KPNO' - Keat Peak National obs. [-111.60 , 31.980]/RAD.
%            'Keck' - Keck observatory.       [-155.478, 19.828]/RAD.
%            'APO'  - Apache Point obs.       [-105.82 , 32.780]/RAD.
%            'CA'   - Calar-Alto obs.         [   2.546, 37.224]/RAD.
%            'MMT'  - MMT obs. (Mt. Hopkins)  [-110.885, 31.688]/RAD.
%            'Paran'- ESO Paranal obs.        [ -70.403,-24.625]/RAD.
%            'LaPal'- La Palma Island obs.    [ -18.881, 28.760]/RAD.
%            'SAO'  - Special Astrop. obs.    [  41.442, 43.653]/RAD.
%            'Palom'- Palomar observatory     [-116.863, 33.357]/RAD.
%            'AAO'  - Anglo-Australian obs.   [ 149.067,-31.277]/RAD.
%            'CTIO' - Cerro Tololo Inter-Am.  [ -70.815,-30.165]/RAD.
%            'ESO'  - ESO La Silla obs.       [ -70.730,-29.257]/RAD.
%          - Date [JD] or [D M Y].
%          - Object R.A. [H M S] or [radians]. 
%          - Object Dec. [Sign D M S] or [radians].
%          - Coordinates equinox, default is 2000.0
% Output : - Structure of plotted parameters, including the following fields:
%            .JD   - Julian day.
%            
% Plot   : Object altitude and airmass during the night, along with
%            the excess in sky brightness, in the object position,
%            due to the Moon illumination.
% Note   : Historically called: obj_obs_cond.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: telescope.obs.daily_observability('Wise',[27 8 2001],[18 0 0],[1 67 0 0]);
% Reliable: 2
%------------------------------------------------------------------------------
RAD      = 180./pi;

SizeRA   = size(RA);
SizeDec  = size(Dec);
SizeDate = size(Date);

if (ischar(ObsCoo)==1)
   [ObsCoo] = telescope.obs.observatory_coo(ObsCoo);
   ObsLon = ObsCoo(1);
   ObsLat = ObsCoo(2);
else
   ObsLon = ObsCoo(:,1);
   ObsLat = ObsCoo(:,2);
end

if (SizeRA(2)==1)
   % do nothing
else
   RA = celestial.coo.convertdms(RA,'H','r');
end
if (SizeDec(2)==1)
   % do nothing
else
   Dec = celestial.coo.convertdms(Dec,'D','R');
end
if (SizeDate(2)==1)
   % already in JD
   JDmid = floor(Date) + 0.5 - ObsLon./(2.*pi);
else
   JDmid = celestial.time.julday([Date, 1-ObsLon(:,1)./(2.*pi)]);
end


JD      = [JDmid-0.5:5./1440:JDmid+0.5].';
UT      = JD+0.5 - floor(JD);
UT      = UT - floor(UT);

N = 12;
JDtick  = [JDmid-0.5:1./N:JDmid+0.5].';
UTtick  = JDtick+0.5 - floor(JDtick);
UTtick  = (UTtick - floor(UTtick)).*24;
UTtickL = ['00000'];
UTtickL = repmat(UTtickL,N,1);
for I=1:1:N
   H = floor(UTtick(I));
   M = round((UTtick(I) - H).*60);
   UTtickL(I,1:5) = sprintf('%2d:%02d',[H,M]);
end

if (nargin==4)
   Equinox = 2000.0;
elseif (nargin==5)
   % do nothing
   if (Equinox~=2000)
      % precess coordinates to J2000.0
      EqStr  = sprintf('j%6.1f',Equinox);
      NewCoo = celestial.coo.coco([RA, Dec],EqStr,'j2000.0','r','r');
      RA     = NewCoo(:,1);
      Dec    = NewCoo(:,2);
   else
      % already in J2000.0
   end
else
   error('Illegal number of input arguments');
end


% calculate sun coordinates:
[SunRA, SunDec]   = celestial.SolarSys.suncoo(JD,'j');
SunHC             = celestial.coo.horiz_coo([SunRA, SunDec],JD,[ObsLon, ObsLat],'h');
SunAz             = SunHC(:,1);
SunAlt            = SunHC(:,2);

% moon coordinates during 24 hours.
[MoonRA, MoonDec] = celestial.SolarSys.mooncool(JD, [ObsLon, ObsLat]);
MoonHC            = celestial.coo.horiz_coo([MoonRA, MoonDec],JD,[ObsLon, ObsLat],'h');
MoonAz            = MoonHC(:,1);
MoonAlt           = MoonHC(:,2);

% Object coordinates during 24 hours
ObjHC             = celestial.coo.horiz_coo([RA, Dec],JD,[ObsLon, ObsLat],'h'); 
ObjAz             = ObjHC(:,1);
ObjAlt            = ObjHC(:,2);
ObjAirMass        = celestial.coo.hardie([pi./2-ObjAlt]);
I = find(ObjAlt<3./RAD);
ObjAirMass(I) = NaN;

% Moon-Object information
C_Ext      = 0.3;
Vsky       = 21.7;
[DelV, MoonElon, ObjMoonDist, MoonIllFrac] = celestial.SolarSys.moon_sky_brightness(JD, [RA, Dec], [ObsLon, ObsLat, 0], C_Ext, Vsky);


% Local Apparent Sidereal Time
LST = celestial.time.lst(JD,ObsLon,'a');

% Sun Rise/Set/Twillight
I_SR = find(diff(sign(SunAlt+0.5./RAD))>0);  % Sun Rise Index
I_SS = find(diff(sign(SunAlt+0.5./RAD))<0);  % Sun Set Index
I_SRT12 = find(diff(sign(SunAlt+12./RAD))>0);  % Sun Morning Twiilight 12 Index
I_SST12 = find(diff(sign(SunAlt+12./RAD))<0);  % Sun Evning Twillight 12 Index
I_SRT18 = find(diff(sign(SunAlt+18./RAD))>0);  % Sun Morning Twiilight 18 Index
I_SST18 = find(diff(sign(SunAlt+18./RAD))<0);  % Sun Evning Twillight 18 Index

% Moon Rise/Set
I_MR = find(diff(sign(MoonAlt+0.5./RAD))>0);  % Moon Rise Index
I_MS = find(diff(sign(MoonAlt+0.5./RAD))<0);  % Moon Set Index

%------------
%--- Plot ---
%------------
plot(JD,ObjAirMass,'k');
%invy;
hold on;
set(gca,'FontSize',12);
h=xlabel('U.T.');
set(h,'FontSize',18);
h=ylabel('Airmass                              ');
set(h,'FontSize',18);

h=text(min(JD)-0.01-0.06,3.0,'/  Moonlight sky-mag excess');
set(h,'rotation',90,'Color',[0 0 1],'FontSize',14);

set(gca,'XTick',JDtick','XTickLabel',UTtickL);
MinY = 0;
MaxY = 5;
axis([min(JD),max(JD),MinY,MaxY]);
% Sun Rise/Set/Twillight
if (isempty(I_SR)==1)
   % do nothing
else
   plot([JD(I_SR);JD(I_SR)],[MinY;MaxY],'r--')
end
if (isempty(I_SS)==1)
   % do nothing
else
   plot([JD(I_SS);JD(I_SS)],[MinY;MaxY],'r--')
end
if (isempty(I_SRT12)==1)
   % do nothing
else
   plot([JD(I_SRT12);JD(I_SRT12)],[MinY;MaxY],'r-.')
end
if (isempty(I_SST12)==1)
   % do nothing
else
   plot([JD(I_SST12);JD(I_SST12)],[MinY;MaxY],'r-.')
end
if (isempty(I_SRT18)==1)
   % do nothing
else
   plot([JD(I_SRT18);JD(I_SRT18)],[MinY;MaxY],'r:')
end
if (isempty(I_SST18)==1)
   % do nothing
else
   plot([JD(I_SST18);JD(I_SST18)],[MinY;MaxY],'r:')
end
% Moon Rise/Set
if (isempty(I_MR)==1)
   % do nothing
else
   plot([JD(I_MR);JD(I_MR)],[MinY;MaxY],'g--')
   %arrow([JD(I_MR) 0.5], [JD(I_MR) 0.0],0.05,0.03);
   hold on;
end
if (isempty(I_MS)==1)
   % do nothing
else
   plot([JD(I_MS);JD(I_MS)],[MinY;MaxY],'g--')
   %arrow([JD(I_MS) 0.0], [JD(I_MS) 0.5],0.05,0.03);
   hold on;
end


% Moon-Object info
plot(JD,abs(DelV),'b');

MoonAM = celestial.coo.hardie([pi./2-MoonAlt]);
I = find(MoonAlt<5./RAD);
MoonAM(I) = NaN;
plot(JD,MoonAM,'b:');


Az = [0:1:90].'./RAD;

AMV = celestial.coo.hardie([pi./2-Az]);

for AMTick=1:0.5:MaxY
   [Min, MinI] = min(abs(AMV-AMTick));
   AzLabel = sprintf('%2d',round(Az(MinI).*RAD));

   h=text(max(JD)+0.01,AMTick,AzLabel);
   set(h,'FontSize',12);
end
h=text(max(JD)+0.01+0.07,3.4,'Altitude, [^{o}]');
set(h,'Rotation',90,'FontSize',18);

Imid = round(length(MoonIllFrac).*0.5);
MoonFrac=sprintf('%+4.2f',mean(MoonIllFrac).*sign(MoonIllFrac(Imid+1)-MoonIllFrac(Imid)));
title(['Moon : ',MoonFrac]);


Res.JD            = JD;
Res.ObjAirMass    = ObjAirMass;
Res.ObjAlt        = ObjAlt;
Res.ObjAz         = ObjAz;
Res.MoonElon      = MoonElon;
Res.ObjMoonDist   = ObjMoonDist;
Res.MoonIllFrac   = MoonIllFrac;
Res.SunAz         = SunAz;
Res.SunAlt        = SunAlt;
Res.MoonAz        = MoonAz;
Res.MoonAlt       = MoonAlt;
Res.MoonExcessSky = DelV;
