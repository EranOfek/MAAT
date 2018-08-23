function Flag=cone_in_polysphere(PolesLong,PolesLat,Long,Lat,Radius)
% Check if a cone (small circle) is within a convex spherical polygon
% Package: celestial.htm
% Description: Check if a cone (small circle) is within a convex spherical
%              polygon which sides are great circles.
% Input  : - Matrix in which each column represent the longitude of the
%            poles of the half-spaces of a spherical polygin, where the
%            pole is directed into the polygon center of mass [rad].
%          - Matrix in which each column represent the latitude of the
%            poles of the half-spaces of a spherical polygin, where the
%            pole is directed into the polygon center of mass [rad].
%          - Vector of longitudes of the cones center [rad].
%            The size is either 1 or like the number of columns of the
%            first and second input arguments.
%          - Vector of latitudes of the cones center [rad].
%            The size is either 1 or like the number of columns of the
%            first and second input arguments.
%          - Vector of radii of the cones [rad].
%            The size is either 1 or like the number of columns of the
%            first and second input arguments.
% Output : - Flag of logical indicating if cone is in polygon.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: HTM=celestial.htm.htm_build(4);
%          Flag=celestial.htm.cone_in_polysphere(HTM(end).PolesCoo(:,1),HTM(end).PolesCoo(:,2),5.5,-0.6,0.01);
%          PLong=rand(3,1000); PLat=rand(3,1000); Long=rand(1,1000); Lat=rand(1000,1); Radius=0.01.*rand(1000,1);
%          Flag=celestial.htm.cone_in_polysphere(PLong,PLat,Long,Lat,Radius);
% Reliable: 2



Long   = Long(:).';
Lat    = Lat(:).';
Radius = Radius(:).';

Dist = acos(bsxfun(@times,sin(PolesLat),sin(Lat)) + bsxfun(@times,cos(PolesLat),cos(Lat)).*cos(bsxfun(@minus,PolesLong,Long)));
Flag = all(bsxfun(@lt,Dist,0.5.*pi+Radius),1);


