function I=inside_celestial_box(RA,Dec,BoxRA,BoxDec,Width,Height)
% Check if coorduinates are within box
% Package: celestial.coo
% Description: Given a list of celestial coordinates, and a box center,
%              width and height, where the box sides are parallel to the
%              coorinate system (i.e., small circles), 
%              check if the coordinates are within the box.
% Input  : - Vector of celestial longitudes [radians] to test.
%          - Vector of celestial latitudes [radians] to test.
%          - Vector of longitudes of center of boxes [radians].
%          - Vector of latitudes of center of boxes [radians].
%          - Vector of full widths of boxes [radians].
%          - Vector of full heights of boxes [radians].
% Output : - Flag indicating if coordinate is inside corresponding box
%            (1) or not (0).
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Mar 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: I=celestial.coo.inside_celestial_box(1,1,1,1,0.1,0.1)
% Reliable: 2
%------------------------------------------------------------------------------

[OffsetLong,OffsetLat,~,~]=celestial.coo.sphere_offset(BoxRA,BoxDec,RA,Dec,'rad','dr');

I = (abs(OffsetLat)<=(0.5.*Height) & abs(OffsetLong)<=(0.5.*Width));

