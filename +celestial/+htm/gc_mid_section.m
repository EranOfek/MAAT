function C=gc_mid_section(Pos1,Pos2,Dim,Output)
% Mid point on great circle between two points
% Package: celestial.htm
% Description: Given two points on a sphere, find the central point lying
%              on the the shortest great circle section connecting the
%              two points.
% Input  : - A list of the first points. This is a three columns matrix
%            in which each row is a unit vector (R1).
%            Alternatively, by setting Dim=2, this can be a three rows
%            matrix in which each column is a unit vector (R1).
%            If R is a two column matrix, then assume the columns are
%            [Long1, Lat1] in radians.
%            If Dim=2, then the input should be a two rows matrix with
%            [Long1, Lat1] in each column.
%          - A list of the second points (similar the the first point).
%          - Dimension along to operate.
%            If 1, then assume the position vectors in (R) and (C) are
%            in rows, while if 2, then assume they are in columns.
%            Default is 1.
%          - Output type: either 'r' for [Long, Lat] in radians,
%            or 'v' for unit vectors. Default is 'v',
% Output : - The coordinates, either unit vectors or [Long, Lat]
%            in radians, of the mid point on the great circle
%            connecting the two points.
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jul 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: C=gc_mid_section([1 1 1],[0 1 0])
% Reliable: 2
%--------------------------------------------------------------------------


Def.Dim    = 1;
Def.Output = 'v';
if (nargin==2)
    Dim = Def.Dim;
    Output = Def.Output;
elseif (nargin==3)
    Output = Def.Output;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end

if (Dim==2)
   Pos1 = Pos1.';
   Pos2 = Pos2.';
end

[P1,Q1]= size(Pos1);
[P2,Q2]= size(Pos2);

if (Q1==2)
    % convert [Long, Lat] to cosine directions
    [CD1, CD2, CD3] = celestial.coo.coo2cosined(Pos1(:,1),Pos1(:,2));
    Pos1 = [CD1, CD2, CD3];
end

if (Q2==2)
    % convert [Long, Lat] to cosine directions
    [CD1, CD2, CD3] = celestial.coo.coo2cosined(Pos2(:,1),Pos2(:,2));
    Pos2 = [CD1, CD2, CD3];
end

% making sure that R and N have the same size
if (P1==1)
    Pos1 = repmat(Pos1,P2,1);
end

if (P2==1)
    Pos2 = repmat(Pos2,P1,1);
end

% find central point by adding vectors
C = Pos1 + Pos2;

% normalize to unit vectors
C = bsxfun(@rdivide,C,sqrt(sum(C.^2,2)));

switch lower(Output)
 case 'r'
    [Long, Lat] = celestial.coo.cosined2coo(C(:,1),C(:,2),C(:,3));
    C = [Long, Lat];
 otherwise
    % do nothing
end

if (Dim==2)
   C = C.';
end
