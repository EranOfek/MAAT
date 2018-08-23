function [RA,Dec,HP]=mooncool(Date,EarthPos,Algo)
% Low-accuracy topocentric equatorial coordinates of the Moon
% Package: celestial.SolarSys
% Description: Calculate low-accuracy topocentric equatorial coordinates
%              of the Moon, referred to the equinox of date.
% Input  : - matrix od dates, [D M Y frac_day] per line,
%            or JD per line. In TT time scale.
%          - [East_Long, North_Lat] of observer in radians.
%            If NaN or empty then calculate geocentric position.
%          - Algorithm:
%            'l' : very low accuracy (default).
%                  0.3 deg in pos. (apparent coordinates).
%                  0.003 deg. in horizontal parallax.
%            'b' : low accuracy ~1' in position.
% Output : - vector of RA, in radians.
%          - vector of Dec. in radians.
%          - Vector of horizontal parallax.
%            r = 1/sin(HP)  SD = 0.2725.*HP
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [RA,Dec,HP]=celestial.SolarSys.mooncool(2451545+(0:1:10)',[1 1]);
% Reliable: 1
%--------------------------------------------------------------------------
RAD = 180./pi;

%FunTPI = inline('(X./(2.*pi) - floor(X./(2.*pi))).*2.*pi','X');


if (nargin==2)
   Algo = 'l';
elseif (nargin==3)
   % do nothing
else
   error('Illigal number of input arguments');
end


SizeDate = size(Date);
N        = SizeDate(1);
ColN     = SizeDate(2);

if (ColN==4)
   JD = celestial.time.julday(Date).';
elseif (ColN==1)
   JD = Date;
else
   error('Illigal number of columns in date matrix');
end

T   = (JD - 2451545.0)./36525.0;

switch Algo
 case 'l'

    n  = JD - 2451545.0;
    if (max(n)>50.*365 || min(n)<-50.*365)
       error('This formulae give good results only between 1950-2050');
    end

    SA1 = sin((134.9 + 477198.85.*T)./RAD);
    SA2 = sin((259.2 - 413335.38.*T)./RAD);
    SA3 = sin((235.7 + 890534.23.*T)./RAD);
    SA4 = sin((269.9 + 954397.70.*T)./RAD);
    SA5 = sin((357.5 +  35999.05.*T)./RAD);
    SA6 = sin((186.6 + 966404.05.*T)./RAD);
    Lam = (218.32 + 481267.883.*T + 6.29.*SA1 - 1.27.*SA2 + 0.66.*SA3 + 0.21.*SA4 - 0.19.*SA5 - 0.11.*SA6)./RAD;
    
    BA1 = sin((93.3  + 483202.03.*T)./RAD);
    BA2 = sin((228.2 + 960400.87.*T)./RAD);
    BA3 = sin((318.3 +   6003.18.*T)./RAD);
    BA4 = sin((217.6 - 407332.20.*T)./RAD);
    Bet = (5.13.*BA1 + 0.28.*BA2 - 0.28.*BA3 - 0.17.*BA4)./RAD;
    
    CA1 = cos((134.9 + 477198.85.*T)./RAD);
    CA2 = cos((259.2 - 413335.38.*T)./RAD);
    CA3 = cos((235.7 + 890534.23.*T)./RAD);
    CA4 = cos((269.9 + 954397.70.*T)./RAD);
    HP  = (0.9508 + 0.0518.*CA1 + 0.0095.*CA2 + 0.0078.*CA3 + 0.0028.*CA4)./RAD;    
    r = 1./sin(HP);
    
    l = cos(Bet).*cos(Lam);
    m = 0.9175.*cos(Bet).*sin(Lam) - 0.3978.*sin(Bet);
    n = 0.3978.*cos(Bet).*sin(Lam) + 0.9175.*sin(Bet);
    
    x = r.*l;
    y = r.*m;
    z = r.*n;

 case 'b'

    [Lon,Lat,Rad,HP] = celestial.SolarSys.moonecool(JD);

    R = 1./sin(HP);
    
    Obl = celestial.coo.obliquity(JD);
    Nt = length(JD);
    x  = zeros(Nt,1);
    y  = zeros(Nt,1);
    z  = zeros(Nt,1);
    for I=1:1:Nt

       L = cos(Lat(I)).*cos(Lon(I));
       M = cos(Obl(I)).*cos(Lat(I)).*sin(Lon(I)) - sin(Obl(I)).*sin(Lat(I));
       N = sin(Obl(I)).*cos(Lat(I)).*sin(Lon(I)) + cos(Obl(I)).*sin(Lat(I));
    
       x(I) = R(I).*L;
       y(I) = R(I).*M;
       z(I) = R(I).*N;   
    end

 otherwise
    error('Unknown algorithm');
end


if (any(isnan(EarthPos)) || isempty(EarthPos))
   % geocentric coordinates
   xt = x;
   yt = y;
   zt = z;
else
   LST = celestial.time.lst(JD,EarthPos(1),'m').*2.*pi;
   Lat = EarthPos(2);

   xt = x - cos(Lat).*cos(LST);
   yt = y - cos(Lat).*sin(LST);
   zt = z - sin(Lat);
end

RA  = atan2(yt,xt);
Dec = asin(zt./sqrt(xt.^2 + yt.^2 + zt.^2));
