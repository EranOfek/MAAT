function [OutPolyX,OutPolyY,SI]=polysort(PolyX,PolyY)
% Sort the vertices of convex polygon by position angle.
% Package: Util.Geom
% Description: Given an (unsorted) convex polygon vertices, sort the
%              vertices acording to their position angle as measured from
%              the polygon center of mass.
% Input  : - Column vector of the X coordinates of the vertices of a
%            polygon.
%          - Column vector of the Y coordinates of the vertices of a
%            polygon.
% Output : - Column vector of the X coordinates of the vertices of a
%            convex polygon.
%          - Column vector of the Y coordinates of the vertices of a
%            convex polygon.
%          - Indices of the input polygon, such that PolyX(I) is the 
%            output polygon.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Mar 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

N = length(PolyX);   %  number of vertices in intersecting polygon

% polygon center of mass
CenterPolyX = mean(PolyX);
CenterPolyY = mean(PolyY);

Theta = atan2(PolyY-CenterPolyY,PolyX-CenterPolyX);
[ST,SI]  = sort(Theta);
OutPolyX = PolyX(SI);
OutPolyY = PolyY(SI);

