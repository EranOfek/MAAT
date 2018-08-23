function [TRS,TAlt,RS_Az]=rise_set(ObjPos,GAST,DefAlt,ObsCoo,DeltaT)
% Calculate rise/set times
% Package: celestial.SolarSys
% Description: Given an object coordinates and observer position,
%              calculate rise/set/transit times, and azimuth and altitude.
%              The times are in the UT1 (not UTC) system.
% Input  : - Object position [RA, Dec] in radians.
%            If object is changing is position, then give:
%            [RA(d-1), Dec(d-1); RA(d), Dec(d); RA(d+1), Dec(d+1)]
%            where RA(x) and Dec(x) are the RA/Dec in:
%            (d-1) day-1 at 0 TDT,
%            (d) day at 0 TDT,
%            (d+1) day+1 at 0 TDT.
%          - Greenwich apparent sidereal time at 0 UT on day (d),
%            in fraction of days.
%          - The geometric altitude of the needed event (radians).
%            e.g., -0.5667/RAD for stars and planets (default).
%                  -0.8333/RAD for the Sun
%                   0.7275*HorizPar - 0.5667/RAD for the Moon.
%                  -6,-12,-18 for twilight.
%          - Observer geodetic position [East Long, North Lat],
%            in radians.
%            Default is the Wise obs. position [  34.763, 30.596]/RAD.
%          - Delta T (=UT-UTC) [fraction of day], default is 0.
% Output : - Time of [Transit; Rise; Set] in fraction of day, for day (d).
%            NaN values for no event. The time is in U.T.
%          - Transit altitude in radians.
%          - Rise/Set azimuth in radians.
% Reference: Meeus J., 1991 in Astronomical Algorithms
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Sep 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%----------------------------------------------------------------------------
RAD = 180./pi;

if (nargin==2),
   DefAlt = -0.5667./RAD;
   ObsCoo = [  34.763, 30.596]./RAD;
   %TimeZone = 2;
   DeltaT = 0;
elseif (nargin==3),
   ObsCoo = [  34.763, 30.596]./RAD;
   %TimeZone = 2;
   DeltaT = 0;
elseif (nargin==4),
   %TimeZone = 2;
   DeltaT = 0;
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

SizePos = size(ObjPos);
if (SizePos(2)~=2),
   error('Object position matrix should have two columns');
end
if (SizePos(1)==1),
   ObjPos = [ObjPos; ObjPos; ObjPos];
elseif (SizePos(1)==3),
   % do nothing
else
   error('Object position matrix should have 1 or 3 raws');
end


DiffRA = diff(ObjPos(:,1));
if (DiffRA(1)<-pi),
   ObjPos([2;3],1) = ObjPos([2;3],1) + 2.*pi;
elseif (DiffRA(2)<-pi),
   ObjPos([3],1) = ObjPos([3],1) + 2.*pi;
elseif (DiffRA(1)>pi),
   ObjPos([2;3],1) = ObjPos([2;3],1) - 2.*pi;
elseif (DiffRA(2)>pi),
   ObjPos([3],1) = ObjPos([3],1) - 2.*pi;
end

%--- patch by Luca Boschini ---
if (abs(ObjPos(1,1))>pi),
   ObjPos(1,1) = ObjPos(1,1) + sign(ObjPos(2,1)).*2.*pi;
else
   if (abs(ObjPos(2,1)-ObjPos(3,1))>pi),
      ObjPos(3,1) = ObjPos(3,1) + sign(ObjPos(2,1)).*2.*pi;
   end
end



if (length(ObsCoo)==2),
   % set Height to 0 meters.
   ObsCoo = [ObsCoo, 0];
end

ObsLon = ObsCoo(1);
ObsLat = ObsCoo(2);

RA  = ObjPos(:,1);
Dec = ObjPos(:,2);

CosH0 = (sin(DefAlt) - sin(ObsLat).*sin(Dec(2)))./(cos(ObsLat).*cos(Dec(2)));
H0  = acos(CosH0);

M    = zeros(3,1);
M(1) = (RA(2) - ObsLon - GAST.*2.*pi)./(2.*pi);    % for the transit [fraction of day].
M(2) = M(1) - H0./(2.*pi);                         % for rising
M(3) = M(1) + H0./(2.*pi);                         % for setting

M    = M - floor(M);

DelM     = ones(3,1);
Accuracy = 0.2./1440;
while (max(abs(DelM))>(Accuracy)),

   % sidereal time at Greenwich [fraction of days]:
   GAST_T = GAST + 360.985647.*M./360;
   GAST_T = GAST_T - floor(GAST_T);

   % interpolate RA/Dec:
   N      = M + DeltaT;

   N1 = N;
   K  = find(abs(N)>1);
   if (isempty(K)==1),
      % do nothing
   else
      N1(K) = sign(N(K));
sprintf('\nWarning - Extrapolating coordinates: %7.4f %7.4f %7.4f\n',N(1),N(2),N(3))
   end
   RA_T   = interp1([-1; 0; 1],RA, N1);
   Dec_T  = interp1([-1; 0; 1],Dec,N1);

   % local HA of object [radians]
   H = GAST_T.*2.*pi + ObsLon - RA_T;
   H = H + pi;
   H = (H./(2.*pi) - floor(H./(2.*pi))).*2.*pi;
   H = H - pi;

   % object altitude and azimuth [radians]
   Alt   = asin(sin(ObsLat).*sin(Dec_T) + cos(ObsLat).*cos(Dec_T).*cos(H));
   SinAz = -cos(Dec_T).*sin(H)./cos(Alt);
   CosAz =  (sin(Dec_T).*cos(ObsLat) - cos(Dec_T).*cos(H).*sin(ObsLat))./cos(Alt);
   Az    = atan2(SinAz, CosAz);

   % corrections [fraction of day]:
   DelM(1)   = -H(1)./(2.*pi);
   DelM(2) = (Alt(2) - DefAlt)./(2.*pi.*cos(Dec_T(2)).*cos(ObsLat).*sin(H(2)));
   DelM(3) = (Alt(3) - DefAlt)./(2.*pi.*cos(Dec_T(3)).*cos(ObsLat).*sin(H(3)));

   M = M + DelM;

end


if (abs(CosH0)>1),
   % object is below/above the horizon all the time
   M(2) = NaN;
   M(3) = NaN;
end

I = find(Az<0);
Az(I) = 2.*pi + Az(I);

TRS   = M;
TAlt  = Alt(1);
RS_Az = [Az(2);Az(3)];


% debug printing
%convertdms(M,'f','H')










