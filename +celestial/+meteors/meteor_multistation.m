function [Az,Alt,Dist,MetCoo]=meteor_multistation(CooList,H,AzAlt)
% Direction for detection of a meteor observed from another station
% Package: celestial.meteors
% Description: Given a list of observers geodetic coordinates in which
%              the first point is a reference point; the azimuth and
%              altitude in which an observer (located in the reference
%              position) is looking to; and the height (H) of a meteor
%              trail - calculate the azimuth and altitude in which
%              observers in other points should look to, in order to
%              detect the same meteor. The function takes into acount
%              the Earth curvature (first order).
% Input  : - CooList : matrix of geodetic coordinates [Long, Lat, Height]
%            in radians and meters, respectively. One raw per observer.
%            The first raw is considered as a reference point.
%            If only two columns are given, the hieght is set
%            to 0 meters.
%          - Meteor height [meters].
%          - [Az, Alt] in radians in which the reference observer
%            is looking to.
% Output : - Vector of azimuths of meteor, per each coordinate [rad].
%          - Vector of altitudes of meteor, per each coordinate [rad].
%          - Observer-meteor distance [meters], per each coordinate.
%          - Coordinates [Long, Lat] for which meteor is in the
%            zenith [radians].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Dec 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Az,Alt,Dist,Coo]=celestial.meteors.meteor_multistation([34.763 30.596; 35 31]./(180./pi),100000,[1 1]);
% Reliable: 2
%--------------------------------------------------------------------------
DATUM = 'WGS84';

CooSize = size(CooList);
if (CooSize(2)==2),
   CooList = [CooList, zeros(CooSize(1),1)];
elseif (CooSize(2)==3),
   % do nothing
else
   error('CooList should contain 2 or 3 columns');
end

RefLong = CooList(1,1);
RefLat  = CooList(1,2);
RefH    = CooList(1,3);
RefAz   = AzAlt(1);
RefAlt  = AzAlt(2);

% calculate earth Radius at observer position [meters]
Data=refellipsoid(DATUM);
A = Data(1);
F = Data(2);
E = sqrt(2.*F - F.^2);
R = A./sqrt(1 - (E.*sin(RefLat)).^2);   % radius of curvature

% Gamma = ref_obs-meteor-earth_center angle
SinGamma = (R+RefH).*sin(0.5.*pi + RefAlt)./(R+H);
Gamma    = asin(SinGamma);
% Delta = ref_obs-earth_center-meteor angle
Delta    = pi - 0.5.*pi - RefAlt - Gamma;
% ref_obs-meteor distance [meters]
if (abs(Delta)<eps),
   Dist1  = H - RefH;
else
   Dist1  = sin(Delta).*(R+RefH)./SinGamma;
end
% Position for which meteor is in the zenith [MetLong, MetLat]
SinMetLat = sin(RefLat).*cos(Delta) + cos(RefLat).*sin(Delta).*cos(RefAz);
MetLat    = asin(SinMetLat);
DL        = asin(sin(Delta).*sin(RefAz)./cos(MetLat));
% DL is a small angle. Therefore, don't need cos(DL)...
MetLong   = RefLong + DL;

MetCoo    = [MetLong, MetLat];

% convert meteor position to cart. coo
[MetGeocCoo,Vm] = celestial.Earth.geod2geoc([MetLong, MetLat, H], DATUM);

Az      = zeros(CooSize(1),1);
Alt     = zeros(CooSize(1),1);
Az(1)   = RefAz;
Alt(1)  = RefAlt;
Dist    = zeros(CooSize(1),1);
Dist(1) = Dist1;
% do for each point in list
for I=2:1:CooSize(1),
   PLong = CooList(I,1);
   PLat  = CooList(I,2);
   PH    = CooList(I,3);
   
   [GeocCoo,V2] = celestial.Earth.geod2geoc(CooList(I,:), DATUM);
   
   % Zeta = earth_center-observer-meteor angle
   V       = V2.*(V2 - Vm);
   AbsV    = sqrt(sum(V2.^2)).*sqrt(sum((V2-Vm).^2));
   CosZeta = sum(V)./AbsV;
   Alt(I)  = acos(CosZeta) - 0.5.*pi;
   Dist(I) = sqrt(sum((V2 - Vm).^2));
   [D,PA]  = celestial.coo.sphere_dist(PLong, PLat, MetLong, MetLat);
   Az(I)   = PA;
end
   
   
   
   
   
   
   
