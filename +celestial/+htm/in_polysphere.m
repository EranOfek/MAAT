function Flag=in_polysphere(Positions,Corners,Crit)
% Is point inside a convex spherical polygon
% Package: celestial.htm
% Description: Check if a list of positions are found within a convex
%              polygon on the celestial sphere in which its sides are
%              great circles.
%              The polygon should be defined according to the
%              right-hand rule.
% Input  : - List of positions to check if they are inside the convex
%            spherical polygon. Each row corrsponds to a position.
%            This is a matrix of either 2 or 3 columns.
%            If two columns are provided then these are [Long, Lat]
%            in radians. If three columns are given, these are
%            cosine directions.
%          - The verteces of a convex polygon on the celestial sphere
%            in which its sides are great circles.
%            Each row correspond to one vertex.
%            This is a matrix of either 2 or 3 columns.
%            If two columns areprovided then these are [Long, Lat]
%            in radians. If three columns are given, these are
%            cosine directions.
%            The coordinates should be ordered such that the
%            right-hand rule is pointing toward the
%            center of the polygon.
%          - Flag indicating if to use ">" or ">=" in in_halfspace.m.
%            If 1 (default) then use (N dot R > C),
%            If 2 then use (N dot R) >= C.
% Output : - A flag indicating if each position (row in Positions
%            matrix) is inside (true) or on/outside the polygon (false).
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jul 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Corners=[0 0;1 0;1 1;0 1];
%          Positions=[0.5 0.5;2 2; 0 0; eps eps];
%          Flag = celestial.htm.in_polysphere(Positions,Corners);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Crit = 1;
if (nargin==2)
   Crit  = Def.Crit;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end

% sort corners
Mean = mean(Corners);
[~,PA] = celestial.coo.sphere_dist(Mean(1),Mean(2),Corners(:,1),Corners(:,2));
[~,SI] = sort(PA,1,'descend');
Corners = Corners(SI,:);

if (size(Positions,2)==2)
   [CD1, CD2, CD3] = celestial.coo.coo2cosined(Positions(:,1),Positions(:,2));
   Positions       = [CD1, CD2, CD3];
end

if (size(Corners,2)==2)
   [CD1, CD2, CD3] = celestial.coo.coo2cosined(Corners(:,1),Corners(:,2));
   Corners         = [CD1, CD2, CD3];
end

% need to sort verteces
[Long,Lat] = celestial.coo.cosined2coo(Corners(:,1),Corners(:,2),Corners(:,3));
[Long,Lat] = celestial.htm.polysphere_sort(Long,Lat);
[CD1,CD2,CD3] = celestial.coo.coo2cosined(Long,Lat);
Corners = [CD1, CD2, CD3];

Corners = [Corners;Corners(1,:)];
% for each pair of verteces call in_halfspace.m
Nvert = size(Corners,1)-1;
% cross (X) two verteces in order to get the pole of the great
% circle defined by the two verteces
PoleVec = Util.math.cross_fast(Corners(1:end-1,:),Corners(2:end,:));
% normalize polar vector to unity
PoleVec = bsxfun(@times,PoleVec, 1./sqrt(sum(PoleVec.^2,2)));
FlagMat = celestial.htm.in_halfspace(Positions,PoleVec,0,1,Crit).';
Flag    = sum(FlagMat,2)==Nvert;
