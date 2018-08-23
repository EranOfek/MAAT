function F=circ_2d(X,Y,Radius,Pos0,Norm)
%--------------------------------------------------------------------------
% circ_2d function                                                 General
% Description: Calculate a circule in 2-D (i.e., cylinder).
% Input  : - Scalar, vector or matrix of X-coordinates in which to calculate
%            the 2-D circle.
%          - same as the x-ccordinates, but for the y-axis.
%          - Radius of circle.
%          - Center of the circle [X, Y].
%            By default Y=X.
%            Default is [0 0].
%            If empty matrix use default.
%          - A normalization constant of the final matrix.
%            If NaN then donot normalize. Default is 1.
% Output : - Value of the 2-D circle in the x-y grid.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Jul 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: gauss_2d.m, triangle_2d.m, lanczos_2d.m
% Example: [MatX,MatY]=meshgrid([-10:1:10],[-10:1:10]);
%          F=circ_2d(MatX,MatY,5,[0 0],1);
%          surface(F);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Pos0   = [0 0];
Def.Norm   = 1;

if (nargin==3),
   Pos0    = Def.Pos0;
   Norm    = Def.Norm;
elseif (nargin==4),
   Norm    = Def.Norm;
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(Pos0)),
   Pos0    = Def.Pos0;
end
if (isempty(Norm)),
   Norm    = Def.Norm;
end

if (length(Pos0)==1),
   Pos0 = [Pos0, Pos0];
end

X0 = Pos0(1);
Y0 = Pos0(2);

F = zeros(size(X));
MatR = sqrt((X-X0).^2 + (Y-Y0).^2);
I = find(MatR<=Radius);
F(I) = 1;

if (isnan(Norm)),
   % do not normalize
else
   F = Norm.*F./sumnd(F);
end
