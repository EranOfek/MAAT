function [X,Y,MatFlag]=cells_intersect_line(X,Y,MatFlag,LineFun)
% Find cells in 2D grid that intersets a line.
% Package: Util.Geom
% Description: Given a grid defining the positions of cells in a 2-D plane
%              and equation of a line in this plane, find all the cells in
%              the grid that intersects the line.
% Input  : - X coordinate of the cells in grid.
%          - Y coordinate of the cells in grid.
%            The X and Y vectors should be uniformly sampled with
%            identical steps.
%          - Matrix of cells in grid, default is zeros(length(Y),length(X)).
%            If empty matrix then use default.
%          - parameters [a b] of line of the form Y=a*x+b;
% Output : - X vector
%          - Y vector
%          - Matrix of flags indicating if the line crossing the cell (1)
%            or not (0).
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Aug 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%----------------------------------------------------------------------------

Xstep = X(2)-X(1);
Ystep = Y(2)-Y(1);  % should be identical to Xstep


if (isempty(MatFlag))
   MatFlag = zeros(length(Y),length(X));
end
[MatX,MatY] = meshgrid(X,Y);
VecX = MatX(:);
VecY = MatY(:);
[MinDist,IP] = dist_p2line(LineFun,[VecX,VecY]);
MatMinDist   = reshape(MinDist,size(MatFlag));
MatFlag      = MatMinDist<0.5.*Xstep;


