function [Geoc,GeocCart]=geod2geoc(Geod,RefEllips)
% Convert Geodetic coordinates to Geocentric coordinates
% Package: celestial.Earth
% Description: Convert Geodetic coordinates to Geocentric coordinates
%              using specified reference ellipsoid.
% Input  : - Geodetic coordinates. This is three column matrix of the
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
% Output : - Geocentric coordinates matrix of the type:
%            [longitude, latitude, radius]
%            where longitude and latitude are measured in radians
%            and radius measured in meters from the reference ellipsoid center.
%          - Geocentric cartesian coordinates [x, y, z] in meters.
% Reference : Astronomical Almnach
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Geoc,GeocCart]=celestial.Earth.geod2geoc([1 1 100]);
% Reliable: 2
%--------------------------------------------------------------------------
if (nargin==1),
   RefEllips  = 'WGS84';
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end
RAD    = 180./pi;
ARCSEC = 1./(3600.*RAD);       % arcsec in radians
N = length(Geod(:,1));
if (length(Geod(1,:))==2),
   Geod = [Geod, zeros(N,1)];
end

RefEllipsData = celestial.Earth.refellipsoid(RefEllips);
A = RefEllipsData(1);
F = RefEllipsData(2);

% Geodetic position
GeodLon = Geod(:,1);
GeodLat = Geod(:,2);
GeodH   = Geod(:,3);
C       = (cos(GeodLat).^2 + ((1-F).*sin(GeodLat)).^2).^(-0.5);
S       = C.*(1 - F).^2;

X       = (A.*C + GeodH).*cos(GeodLat).*cos(GeodLon);
Y       = (A.*C + GeodH).*cos(GeodLat).*sin(GeodLon);
Z       = (A.*S + GeodH).*sin(GeodLat);

R2 = X.^2 + Y.^2;    % =(GeocRad.*cos(GeoLat)).^2;
GeocLat = atan(Z./sqrt(R2));
GeocRad = sqrt(R2)./cos(GeocLat);

Geoc = [GeodLon, GeocLat, GeocRad];

GeocCart = [X, Y, Z];
    
      
