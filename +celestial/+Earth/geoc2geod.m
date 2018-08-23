function Geod=geoc2geod(Geoc,RefEllips,Accuracy)
% Convert Geocentric coordinates to Geodetic coordinates
% Description: Convert Geocentric coordinates to Geodetic coordinates
%              using specified reference ellipsoid.
% Input  : - Geocentric coordinates. This is three column matrix of the
%            type: [Longitude, Latitude, Radius]
%            where Radius is the distance from the sphere center.
%            Units: Longitude and latitude measured in radians while
%            the radius is measured in meters above reference ellipsoid.
%          - Reference ellipsoid:
%            'Merit1983' - Merit 1983           a=6378137  1/f=298.257
%            'GRS80'     - GRS 80 (IUGG 1980)     6378137      298.257222
%            'GRS67'     - GRS 67 (IUGG 1967)     6378160      298.247167
%            'IAU1976'   - IAU 1976               6378140      298.257
%            'IAU1964'   - IAU 1964               6378160      298.25
%            'WGS84'     - WGS 1984 (default)
%          - Convergence accuracy in arcsec, (default is 0.01 arcsec).
%            (Iterative algorithm).
% Output : - Geodetic coordinates matrix of the type:
%            [longitude, latitude, Height]
%            where longitude and latitude are measured in radians
%            and Height measured in meters above the reference ellipsoid.
% Reference : Astronomical Almnach
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Geod=celestial.Earth.geoc2geod([1 1 7000]);
% Reliable: 2
%------------------------------------------------------------------------------
if (nargin==1),
   RefEllips  = 'WGS84';
   Accuracy   = 0.01;    % arcsec   
elseif (nargin==2),
   Accuracy   = 0.01;    % arcsec   
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(RefEllips)==1),
   RefEllips = 'WGS84';
end

RAD    = 180./pi;
ARCSEC = 1./(3600.*RAD);       % arcsec in radians

N = length(Geoc(:,1));


RefEllipsData = celestial.Earth.refellipsoid(RefEllips);
A = RefEllipsData(1);
F = RefEllipsData(2);

% Geodetic position
GeocLon = Geoc(:,1);
GeocLat = Geoc(:,2);
GeocRad = Geoc(:,3);

% initial guess
GeodLon = GeocLon;

X       = GeocRad.*cos(GeocLat).*cos(GeocLon);
Y       = GeocRad.*cos(GeocLat).*sin(GeocLon);
Z       = GeocRad.*sin(GeocLat);

R       = sqrt(X.^2 + Y.^2);

E2 = 2.*F - F.^2;

Phi1 = atan(Z./R);

DelPhi = 1;

while (max(DelPhi)>(Accuracy.*ARCSEC)),
   C    = (1 - E2.*sin(Phi1).^2).^(-0.5);
   Phi  = atan((Z + A.*C.*E2.*sin(Phi1))./R);
   H    = R./cos(Phi) - A.*C;
   
   DelPhi = abs(Phi-Phi1);
   Phi1   = Phi; 
end

Geod = [GeodLon, Phi, H];
