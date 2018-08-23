function [Flag,Dist,PA]=sphere_dist_thresh(Long,Lat,LongRef,LatRef,MaxDist,Shape)
%--------------------------------------------------------------------------
% sphere_dist_thresh function                                        ephem
% Description: Given Long and Lat coordinates and a reference coordinates
%              (in radians) return a flag indicating if each point is
%              within a spherical distance from a reference point.
% Input  : - Matrix of Longitude coordinates [rad].
%          - Matrix of Latitude coordinates [rad].
%          - Reference Longitude coordinates [rad].
%          - Reference Latitude coordinates [rad].
%          - Distance threshold [rad] (radius, or box full width).
%          - Search shape {'box'|'circ'}. Default is 'circ'.
% Output : - Flag indicating if the input Long/Lat points are in the
%            requested region.
%          - Distance.
%          - Position angle [radians].
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [Flag,Dist,PA]=sphere_dist_thresh(rand(100,1),rand(100,1),0.5,0.5,0.2,'box')
% Reliable: 2
%--------------------------------------------------------------------------
RAD  = 180./pi;
Crit = 2;

if (nargin<6)
    Shape = 'circ';
end

switch lower(Shape)
    case 'circ'
        [Dist,PA] = celestial.coo.sphere_dist(Long,Lat,LongRef,LatRef);
        Flag      = Dist<=MaxDist;
    case 'box'
        [LatCor,LongCor] = reckon(LatRef.*RAD,LongRef.*RAD,sqrt(2).*0.5.*MaxDist.*RAD, 45 +(0:90:270)');
        Corners   = [LongCor, LatCor]./RAD;
        Flag      = celestial.htm.in_polysphere([Long Lat],Corners,Crit);
        
        if (nargout>2)
            [Dist,PA] = celestial.coo.sphere_dist(Long,Lat,LongRef,LatRef);
        end
    otherwise
        error('Unknown Shape option');
end
