function [G,Gt]=topocentric_vec(JD_UT1, Geod, RefEllips, Xp, Yp);
%--------------------------------------------------------------------------
% topocentric_vec function                                           ephem
% Description: Calculate the topocentric position and velocity vectors
%              of an observer, with respect to the true equator and
%              equinox of date. In otder to transform the vectors to a
%              coodinates system in respect to the Earth's mean equator
%              and equinox, use inv(P)*inv(N)*G and inv(P)*inv(N)*Gt,
%              where P and N are the precession and nutation rotation
%              matrices (see rotm_coo.m).
% Input  : - Column vector of JD in UT1 time scale,
%            or date [D M Y Frac] | [D M Y H M S] | [D M Y].
%          - Geodetic coordinates. This is three column matrix of the
%            type: [Longitude, Latitude, Height], in case that two
%            column matrix is given, than the height is taken as zero.
%            Units: Longitude and latitude measured in radians while
%            the height is measured in meters above reference ellipsoid.
%          - Reference ellipsoid:
%            'Merit1983' - Merit 1983           a=6378137  1/f=298.257
%            'GRS80'     - GRS 80 (IUGG 1980)     6378137      298.257222
%            'GRS67'     - GRS 67 (IUGG 1967)     6378160      298.247167
%            'IAU1976'   - IAU 1976               6378140      298.257
%            'IAU1964'   - IAU 1964               6378160      298.25
%            'WGS84'     - WGS 1984 (default)
%          - The X coordinate of the celestial ephemeris pole
%            with respect to the terrestrial pole measured along
%            the 0 meridian (radians). - default is 0.
%          - The Y coordinate of the celestial ephemeris pole
%            with respect to the terrestrial pole measured along
%            the 270 meridian (radians). default is 0.
% Output : - The topocentric position vector with respect to the
%            true equator and equinox of date. [meters]
%            Each column for one time.
%          - The topocentric velocity vector with respect to the
%            true equator and equinox of date. [meters/sec]
%            Each column for one time.
% Tested : Matlab 5.3
%     By : Eran O. Ofek             May 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html 
%--------------------------------------------------------------------------

RADIAN = 180./pi;
N = size(JD_UT1,1);
W = 7.2921151467e-5;  %rad/sec (Earth ang. velocity)

if (nargin==2),
   RefEllips = 'WGS84';
   Xp = zeros(N,1);
   Yp = zeros(N,1);
elseif (nargin==3),
   Xp = zeros(N,1);
   Yp = zeros(N,1);
elseif (nargin==4),
   Yp = zeros(N,1);
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (size(JD_UT1,2)==1),
   % do nothing
else
   % convert date to JD
   JD_UT1 = julday(JD_UT1).';
end


[Geoc,GeocCart]=geod2geoc(Geod,RefEllips);
GeocCart = GeocCart.';
if (size(GeocCart,2)==1),
   GeocCart = GeocCart*ones(1,N);
end

LAST = lst(JD_UT1, 0, 'a');     % calculate app. sidereal time at Greenwich
LAST = LAST.*2.*pi;             % convert to radians

G  = zeros(3,N);
Gt = zeros(3,N);

for I=1:1:N,
   R = GeocCart(:,I);
   if (isnan(Xp(I))==1),
      R2Xp = diag([1 1 1]);
   else
      R2Xp = rotm(Xp(I),2);
   end

   if (isnan(Yp(I))==1),
      R1Yp = diag([1 1 1]);
   else
      R1Yp = rotm(Yp(I),1);
   end

   G(:,I) = rotm(+LAST(I),3)*R1Yp*R2Xp*R;
   RotG   = [-sin(LAST(I)) -cos(LAST(I)) 0; cos(LAST(I)) -sin(LAST(I)) 0; 0 0 0];
   Gt(:,I) = W.*RotG*R1Yp*R2Xp*R;
   %  [-G(2).*W; G(1).*W; 0];
end

