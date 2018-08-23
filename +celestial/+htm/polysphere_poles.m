function [PoleLong,PoleLat]=polysphere_poles(VertLong,VertLat)
% Given a spherical polygon vertces, find the poles of each of sides
% Package: celestial.htm
% Description: Given the longitude and latitude of a convex polygon on the
%              sphere in which its sides are great circles, find the poles
%              of each great circle aiming the polygon center of mass.
% Input  : - Vector of the longitude of the vertces [radians].
%          - Vector of the latitude of the vertces [radians].
% Output : - Vector of poles longitude [radians].
%          - Vector of poles latitude [radians].
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VertLong=[0 1 0]'; VertLat=[0 0 1]';
%          [PoleLong,PoleLat]=celestial.htm.polysphere_poles(VertLong,VertLat)
% Reliable: 2


VertLong = VertLong(:);
VertLat  = VertLat(:);
% calculate the polygon center of mass
[CD1,CD2,CD3]    = celestial.coo.coo2cosined(VertLong,VertLat);
[CenLong,CenLat] = celestial.coo.cosined2coo(mean(CD1),mean(CD2),mean(CD3));
% calcualate direction via P.A.
[~,PA] = celestial.coo.sphere_dist_fast(CenLong,CenLat,VertLong,VertLat);
[~,SI] = sort(PA);
VertLong = VertLong(SI);
VertLat  = VertLat(SI);
CD1 = CD1(SI);
CD2 = CD2(SI);
CD3 = CD3(SI);
Corners = [[CD1, CD2, CD3]; [CD1(1), CD2(1), CD3(1)]];
PoleVec = Util.math.cross_fast(Corners(1:end-1,:),Corners(2:end,:));
% convert to long/lat
[PoleLong,PoleLat] = celestial.coo.cosined2coo(PoleVec(:,1),PoleVec(:,2),PoleVec(:,3));
