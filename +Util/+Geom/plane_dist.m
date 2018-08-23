function [Dist,PA]=plane_dist(X1,Y1,varargin)
% Distance between points on a 2D plane.
% Package: Util.Geom
% Description: Calculate the planner distance and angle between points.
% Input  : - Column vector of X1 coordinates.
%            If only two input arguments are provided than this a matrix of
%            the [X1,Y1] coordinates.
%          - Column vector of Y1 coordinates.
%            If only two input arguments are provided than this a matrix of
%            the [X2,Y2] coordinates.
%          - Column vector of X2 coordinates.
%          - Column vector of Y2 coordinates.
% Output : - Column vector of distance between pairs of points.
%          - Angles between pairs of points.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: sphere_dist.m
% Example: [Dist,PA]=Util.Geom.plane_dist(1,1,2,2);
%          [Dist,PA]=Util.Geom.plane_dist([1 1;2 1],[1 2; 2 2]);
% Reliable: 1
%--------------------------------------------------------------------------
if (nargin==2)
   Y2 = Y1(:,2);
   X2 = Y1(:,1);

   Y1 = X1(:,2);
   X1 = X1(:,1);
else
   X2 = varargin{1};
   Y2 = varargin{2};
end

Dist = sqrt((X1-X2).^2 + (Y1-Y2).^2);
if (nargout>1)
   PA   = atan2((Y2-Y1),(X2-X1));
end
