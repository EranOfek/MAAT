function [MinDist,IP]=dist_p2line(Line,P0)
% Distance between point and a line.
% Package: Util.Geom
% Description: Calculate the minimum distance in a 2-d space between a line
%              and a point.
% Input  : - Line definition,
%            If a two element vector is given [a b] then the line is
%            defined by y=ax+b;
%            If 2x2 matrix [x1 x2; y1 y2] is given then the line is
%            defined as going through the points (x1,y1) and (x2,y2).
%          - Point position [xp,yp].
% Output : - Minimum distance between line and point.
%          - The coordinates (x,y) of the nearest-point-on-the-line
%            to the point [xp,yp].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Feb 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [MinDist,IP]=Util.Geom.dist_p2line([1,0],[1,0]);
% Reliable: 2
%--------------------------------------------------------------------------
[Ni,Nj] = size(Line);

if (Ni.*Nj==2)
   P1 = [0; Line(1).*0 + Line(2)];
   P2 = [1; Line(1).*1 + Line(2)];
elseif (Ni.*Nj==4)
   P1 = Line(:,1);
   P2 = Line(:,2);
else
   error('Illegal Line definition');
end

Np = size(P0,1);
IP = size(Np,2);
MinDist = size(Np,1);
for I=1:1:Np
   U  = ((P0(I,1)-P1(1)).*(P2(1)-P1(1)) + (P0(I,2)-P1(2)).*(P2(2)-P1(2)))./ ...
        ((P1(1)-P2(1)).^2 + (P1(2)-P2(2)).^2);

   % intersection point
   IP(I,1:2) = [P1(1) + U.*(P2(1)-P1(1)),    P1(2) + U.*(P2(2)-P1(2))];

   MinDist(I) = sqrt((P0(I,1)-IP(I,1)).^2 + (P0(I,2)-IP(I,2)).^2);
end
