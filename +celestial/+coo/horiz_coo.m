function [OutCoo,Lat]=horiz_coo(InCoo,JD,TopoPos,Direction)
% Celestial equatorial coordinates to horizontal coordinates
% Package: celestial.coo
% Description: Convert Right Ascension and Declination to horizontal
%              coordinates or visa versa.
% Input  : - Two columns matrix of coordinates,
%            (Long & Lat) | (Az & Alt) in radians
%          - vector of JDs + UT fraction of day,
%            if scalar value is given, then it duplicated
%            for all coordinates.
%          - Geodetic Coordinates, east long & north lat in radians
%            if scalar value is given then it is duplicate for
%            all coordinates.
%          - Direction,
%            'h' - from equatorial to horizontal (default).
%            'e' - from horizontal to equatorial.
% Output : - two column matrix of output coordinates.
%          * If two output arguments are requested than return the two
%            coordinates seperatly.
% Tested : matlab 5.3
%     By : Eran O. Ofek                    Aug 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: OutCoo=celestial.coo.horiz_coo([1 1],celestial.time.julday([1 1 2015]),[1 1])
% Reliable: 1
%--------------------------------------------------------------------------
if (nargin==3)
   Direction = 'h';
elseif (nargin==4)
   % no default
else
   error('Illigal number of input arguments');
end

N = length(InCoo(:,1));
if (length(JD)==1)
   JD = JD.*ones(N,1);
end

if (numel(TopoPos(:,1))==1)
   TopoPos = ones(N,1)*TopoPos;
end


% Don't convert Geodetic latitude to Geocentric.
%GeodLat = TopoPos(:,2);
%GeocLatTemp = geod2geoc([TopoPos(:,1),GeodLat,zeros(N,1)],'WGS84');
%GeocLat     = GeocLatTemp(:,2);
%TopoPos(:,2) = GeocLat;


% calculating Local Mean Sidereal Time
LST=celestial.time.lst(JD,TopoPos(:,1),'m');


if (Direction=='h')
   % convert equatorial to horizontal

   % calculate the Hour Angle
   HA  = LST.*2.*pi - InCoo(:,1);

   Dec = InCoo(:,2);
   Lat = TopoPos(:,2);

   SinAlt = sin(Dec).*sin(Lat) + cos(Dec).*cos(HA).*cos(Lat);
   CosAlt = sqrt(1-SinAlt.*SinAlt);

   SinAz  = (-cos(Dec).*sin(HA))./CosAlt;
   CosAz  = (sin(Dec).*cos(Lat) - cos(Dec).*cos(HA).*sin(Lat))./CosAlt;

   Az     = atan2(SinAz, CosAz);
   Alt    = asin(SinAlt);

   I      = find(Az<0);
   Az(I)  = 2.*pi+Az(I);

   OutCoo = [Az, Alt];
elseif (Direction=='e')
   % horizontal to equatorial
   Az     = InCoo(:,1);
   Alt    = InCoo(:,2);
   Lat    = TopoPos(:,2);

   SinDec = sin(Alt).*sin(Lat) + cos(Alt).*cos(Az).*cos(Lat);
   CosDec = sqrt(1 - SinDec.*SinDec);

   SinHA  = (-cos(Alt).*sin(Az))./CosDec;
   CosHA  = (sin(Alt).*cos(Lat) - cos(Alt).*cos(Az).*sin(Lat))./CosDec;
   HA     = atan2(SinHA, CosHA);

   RA     = LST.*2.*pi - HA;
   Dec    = asin(SinDec);

   % converting to range [0,1)
   RA     = 2.*pi.*(RA./(2.*pi) - floor(RA./(2.*pi)));

   OutCoo = [RA, Dec];
else
   error('Illigal type of conversion');
end



if nargout>1
    Lat = OutCoo(:,2);
    OutCoo  = OutCoo(:,1);
end







