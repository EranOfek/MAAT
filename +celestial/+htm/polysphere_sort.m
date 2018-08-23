function [OutPolyLong,OutPolyLat]=polysphere_sort(PolyLong,PolyLat)
% Sort a convex spherical polygon
% Package: celestial.htm
% Description: Given an (unsorted) convex spherical polygon vertices, sort the
%              vertices acording to their descending position angle as
%              measured from the polygon center of mass.
% Input  : - Column vector of the Long coordinates of the vertices of a
%            polygon [rad].
%          - Column vector of the Lat coordinates of the vertices of a
%            polygon [rad].
% Output : - Column vector of the Long coordinates of the vertices of a
%            convex polygon [rad].
%          - Column vector of the Lat coordinates of the vertices of a
%            convex polygon [rad].
%          - Indices of the input polygon, such that PolyX(I) is the 
%            output polygon.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Sep 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [OutPolyLong,OutPolyLat]=celestial.htm.polysphere_sort([0 1 2 1]',[0 0 1 1]');
% Reliable: 2
%--------------------------------------------------------------------------

[CD1,CD2,CD3] = celestial.coo.coo2cosined(PolyLong,PolyLat);
[MeanLong, MeanLat] = celestial.coo.cosined2coo(mean(CD1),mean(CD2),mean(CD3));

[~,PA] = celestial.coo.sphere_dist(MeanLong,MeanLat,PolyLong,PolyLat);
[~,SI] = sort(PA,1,'descend');
OutPolyLong = PolyLong(SI);
OutPolyLat = PolyLat(SI);
