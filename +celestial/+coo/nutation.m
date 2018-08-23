function [Nut,NutMatrix]=nutation(JD,MatType)
% Intermidiate accuracy IAU 1984 nutation
% Package: celestial.coo
% Description: Calculate the Nutation in longitude and latitude, and the
%              nutation rotation matrix.
%              This is a low accuracy version based on the IAU 1984 nutation
    %              series. See also: nutation1984.m
% Input  : - Column vector of julian day.
%          - Type of nutation matrix:
%            'f' : full precision, (default).
%            'l' : linearized matrix.
% Output : - Matrix listing nutation in longitude and obliquity,
%            line per JD [radians], in the first and second columns,
%            respectively.
%          - Nutation rotation matrix.
%            In case that more then one JD is asked,
%            then this is a cube.
% Reference : Explanatory Supplement to the Astronomical Almanac (1992)
% See also: nutation1984.m, nutation2rotmat.m
% Tested : Matlab 5.2
%     By : Eran O. Ofek                    Feb 2000
%    URL : htpp://weizmann.ac.il/home/eofek/matlab/
% Example: [N,NM]=celestial.coo.nutation(2451545);
% Reliable: 2
%--------------------------------------------------------------------------
if (nargin==1)
   MatType = 'f';
elseif (nargin==2)
   % no default
else
   error('Illegal number of input arguments');
end

RADIAN = 180./pi;

T   = (JD - 2451545.0)./36525.0;

% Mean elongation of the Moon from the Sun:
D = (297.85036 + 445267.111480.*T - 0.0019142.*T.*T + T.*T.*T./189474)./RADIAN;

% Mean anomaly of the Sun (Earth):
M = (357.52772 + 35999.050340.*T - 0.0001603.*T.*T - T.*T.*T./300000)./RADIAN;

% Mean anomaly of the Moon:
Mt = (134.96298 + 477198.867398.*T + 0.0086972.*T.*T + T.*T.*T./56250)./RADIAN;

% Moon's argument of latitude
F = (93.27191 + 483202.017538.*T - 0.0036825.*T.*T + T.*T.*T./327270)./RADIAN;

% Longitude of ascending node of the Moon's mean orbit
Om = (125.04452 - 1934.136261.*T + 0.0020708.*T.*T + T.*T.*T./450000)./RADIAN;

% nutation in longitude [ 0."0001]
DLon = (-171996 - 174.2.*T).*sin(Om) + ...,
       ( -13187 +   1.6.*T).*sin(2.*(-D + F + Om)) + ...,
       (  -2274 -   0.2.*T).*sin(2.*(F + Om)) + ...,
       (   2062 +   0.2.*T).*sin(2.*Om) + ...,
       (   1426 -   3.4.*T).*sin(M) + ...,
       (    712 +   0.1.*T).*sin(Mt) + ...,
       (   -517 +   1.2.*T).*sin(2.*(-D + F + Om) + M) + ...,
       (   -386 -   0.4.*T).*sin(2.*F + Om) + ...,
       (   -301           ).*sin(Mt + 2.*(F + Om)) + ...,
       (    217 -   0.5.*T).*sin(-M + 2.*(-D + F + Om)) + ...,
       (   -158           ).*sin(-2.*D + Mt) + ...,
       (    129 +   0.1.*T).*sin(2.*(-D + F) + Om) + ...,
       (    123           ).*sin(-Mt + 2.*(F + Om)) + ...,
       (     63           ).*sin(2.*D) + ...,
       (     63 +   0.1.*T).*sin(Mt + Om) + ...,
       (    -59           ).*sin(-Mt + 2.*(D + F + Om)) + ...,
       (    -58 -   0.1.*T).*sin(-Mt + Om) + ...,
       (    -51           ).*sin(Mt + 2.*F + Om) + ...,
       (     48           ).*sin(2.*(-D + Mt)) + ...,
       (     46           ).*sin(Om + 2.*(-Mt + F)) + ...,
       (    -38           ).*sin(2.*(D + F + Om)) + ...,
       (    -31           ).*sin(2.*(Mt + F + Om)) + ...,
       (     29           ).*sin(2.*Mt) + ...,
       (     29           ).*sin(Mt + 2.*(-D + F + Om)) + ...,
       (     26           ).*sin(2.*F) + ...,
       (    -22           ).*sin(2.*(-D + F)) + ...,
       (     21           ).*sin(-Mt + Om + 2.*F) + ...,
       (     17 -   0.1.*T).*sin(2.*M) + ...,
       (     16           ).*sin(2.*D - Mt + Om) + ...,
       (    -16 +   0.1.*T).*sin(2.*(-D + M + F + Om)) + ...,
       (    -15           ).*sin(M + Om) + ...,
       (    -13           ).*sin(-2.*D + Mt + Om) + ...,
       (    -12           ).*sin(-M + Om) + ...,
       (     11           ).*sin(2.*(Mt - F)) + ...,
       (    -10           ).*sin(2.*(D + F) - Mt + Om) + ...,
       (     -8           ).*sin(2.*(D + F + Om) + Mt) + ...,
       (      7           ).*sin(2.*(F + Om) + M) + ...,
       (     -7           ).*sin(-2.*D + M + Mt) + ...,
       (     -7           ).*sin(-M + 2.*(F + Om)) + ...,
       (     -7           ).*sin(2.*(D + F) + Om) + ...,
       (      6           ).*sin(2.*D + Mt) + ...,
       (      6           ).*sin(2.*(-D + Mt + F + Om)) + ...,
       (      6           ).*sin(2.*(-D + F) + Mt + Om) + ...,
       (     -6           ).*sin(2.*(D - Mt) + Om) + ...,
       (     -6           ).*sin(2.*D + Om) + ...,
       (      5           ).*sin(-M + Mt) + ...,
       (     -5           ).*sin(2.*(F - D) + Om - M) + ...,
       (     -5           ).*sin(Om - 2.*D) + ...,
       (     -5           ).*sin(2.*(Mt + F) + Om) + ...,
       (      4           ).*sin(2.*(Mt - D) + Om) + ...,
       (      4           ).*sin(2.*(F - D) + M + Om) + ...,
       (      4           ).*sin(Mt - 2.*F) + ...,
       (     -4           ).*sin(Mt - D) + ...,
       (     -4           ).*sin(M -2.*D) + ...,
       (     -4           ).*sin(D) + ...,
       (      3           ).*sin(Mt + 2.*F) + ...,
       (     -3           ).*sin(2.*(F + Om - Mt)) + ...,
       (     -3           ).*sin(Mt - D - M) + ...,
       (     -3           ).*sin(M + Mt) + ...,
       (     -3           ).*sin(Mt - M + 2.*(F - Om)) + ...,
       (     -3           ).*sin(2.*(D + F + Om) - M - Mt) + ...,
       (     -3           ).*sin(3.*Mt + 2.*(F + Om)) + ...,
       (     -3           ).*sin(2.*(D + F + Om) - M);
       
    
    
% nutation in obliquity [ 0."0001]
DObl = (  92025 +   8.9.*T).*cos(Om) + ...,
       (   5736 -   3.1.*T).*cos(2.*(-D + F + Om)) + ...,
       (    977 -   0.5.*T).*cos(2.*(F + Om)) + ...,
       (   -895 +   0.5.*T).*cos(2.*Om) + ...,
       (     54 -   0.1.*T).*cos(M) + ...,
       (     -7           ).*cos(Mt) + ...,
       (    224 -   0.6.*T).*cos(2.*(-D + F + Om) + M) + ...,
       (    200           ).*cos(2.*F + Om) + ...,
       (    129 -   0.1.*T).*cos(Mt + 2.*(F + Om)) + ...,
       (    -95 +   0.3.*T).*cos(-M + 2.*(-D + F + Om)) + ...,
       (    -70           ).*cos(2.*(-D + F) + Om) + ...,
       (    -53           ).*cos(-Mt + 2.*(F + Om)) + ...,
       (    -33           ).*cos(Mt + Om) + ...,
       (     26           ).*cos(-Mt + 2.*(D + F + Om)) + ...,
       (     32           ).*cos(-Mt + Om) + ...,
       (     27           ).*cos(Mt + 2.*F + Om) + ...,
       (    -24           ).*cos(Om + 2.*(-Mt + F)) + ...,
       (     16           ).*cos(2.*(D + F + Om)) + ...,
       (     13           ).*cos(2.*(Mt + F + Om)) + ...,
       (    -12           ).*cos(Mt + 2.*(-D + F + Om)) + ...,
       (    -10           ).*cos(-Mt + Om + 2.*F) + ...,
       (     -8           ).*cos(2.*D - Mt + Om) + ...,
       (      7           ).*cos(2.*(-D + M + F + Om)) + ...,
       (      9           ).*cos(M + Om) + ...,
       (      7           ).*cos(-2.*D + Mt + Om) + ...,
       (      6           ).*cos(-M + Om) + ...,
       (      5           ).*cos(2.*(D + F) - Mt + Om) + ...,
       (      3           ).*cos(2.*(D + F + Om) + Mt) + ...,
       (     -3           ).*cos(2.*(F + Om) + M) + ...,
       (      3           ).*cos(-M + 2.*(F + Om)) + ...,
       (      3           ).*cos(2.*(D + F) + Om) + ...,
       (     -3           ).*cos(2.*(-D + Mt + F + Om)) + ...,
       (     -3           ).*cos(2.*(-D + F) + Mt + Om) + ...,
       (      3           ).*cos(2.*(D - Mt) + Om) + ...,
       (      3           ).*cos(2.*D + Om) + ...,
       (      3           ).*cos(2.*(F - D) + Om - M) + ...,
       (      3           ).*cos(Om - 2.*D) + ...,
       (      3           ).*cos(2.*(Mt + F) + Om);

    
% convert to radians
DLon = DLon.*0.0001./(3600.*RADIAN);
DObl = DObl.*0.0001./(3600.*RADIAN);
Nut = [DLon, DObl];
    
% calculate nutation rotation matrix
if (nargout>1)
   NutMatrix = ones(3,3,length(JD));
   Obl = celestial.coo.obliquity(JD);
   switch MatType
    case {'f'}
       % Full precesion Nutation Matrix
       Obl1 = Obl + DObl;
       NutMatrix(1,1,:) = cos(DLon);
       NutMatrix(1,2,:) = -sin(DLon).*cos(Obl);
       NutMatrix(1,3,:) = -sin(DLon).*sin(Obl);
       NutMatrix(2,1,:) = sin(DLon).*cos(Obl1);
       NutMatrix(2,2,:) = cos(DLon).*cos(Obl1).*cos(Obl) + sin(Obl1).*sin(Obl);
       NutMatrix(2,3,:) = cos(DLon).*cos(Obl1).*sin(Obl) - sin(Obl1).*cos(Obl);
       NutMatrix(3,1,:) = sin(DLon).*sin(Obl1);
       NutMatrix(3,2,:) = cos(DLon).*sin(Obl1).*cos(Obl) - cos(Obl1).*sin(Obl);
       NutMatrix(3,3,:) = cos(DLon).*sin(Obl1).*sin(Obl) + cos(Obl1).*cos(Obl);
    case {'l'}
       % Linearized Nutation matrix
       %NutMatrix(1,1,:) = 1;
       NutMatrix(1,2,:) = -DLon.*cos(Obl);
       NutMatrix(1,3,:) = -DLon.*sin(Obl);
       NutMatrix(2,1,:) = DLon.*cos(Obl);
       %NutMatrix(2,2,:) = 1;
       NutMatrix(2,3,:) = -DObl; 
       NutMatrix(3,1,:) = DLon.*sin(Obl);
       NutMatrix(3,2,:) = DObl;
       %NutMatrix(3,3,:) = 1;
    otherwise
       error('Unknown MatType');
   end
end
