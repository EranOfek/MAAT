function [Flag,Dist,PA]=plane_dist_thresh(X,Y,Xr,Yr,MaxDist,Shape)
% Check if the distance between points on a plane is below some value.
% Package: Util.Geom
% Description: Given X and Y coordinates and a reference coordinates
%              return a flag indicating if each point is within a distance
%              from a reference point (using planer distance).
% Input  : - Matrix of X coordinates.
%          - Matrix of Y coordinates.
%          - Reference X coordinates.
%          - Reference Y coordinates.
%          - Distance threshold.
%          - Search shape {'box'|'circ'}. Default is 'circ'.
% Output : - Flag indicating if the input X/Y points are in the requested
%            region.
%          - Distance.
%          - Position angle [radians].
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Flag,Dist,PA]=Util.Geom.plane_dist_thresh(rand(100,1),rand(100,1),0.5,0.5,0.2,'box')
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<6)
    Shape = 'circ';
end

Dist = sqrt((X-Xr).^2 + (Y-Yr).^2);

switch lower(Shape)
    case 'box'
        Flag = abs(X-Xr)<=MaxDist & abs(Y-Yr)<=MaxDist;
    case 'circ'
        Flag = Dist<=MaxDist;
    otherwise
        error('Unknown Shape option');
end

if (nargout>2)
   PA   = atan2((Yr-Y),(Xr-X));
end
