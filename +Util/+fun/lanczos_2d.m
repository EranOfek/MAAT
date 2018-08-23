function F=lanczos_2d(X,Y,A,Stretch,Pos0,Norm)
%--------------------------------------------------------------------------
% lanczos_2d function                                              General
% Description: Calculate lanczos function in 2-D.
% Input  : - Scalar, vector or matrix of X-coordinates in which to calculate
%            the 2-D triangle.
%          - same as the x-ccordinates, but for the y-axis.
%          - The order of the lanczos function. Default is 2.
%          - Radial stretch factor. Default is 1.
%          - Center of the lanczos function [X, Y].
%            By default Y=X.
%            Default is [0 0].
%            If empty matrix use default.
%          - A normalization constant of the final matrix.
%            If NaN then donot normalize. Default is 1.
% Output : - Value of the 2-D box in the x-y grid.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Jul 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: gauss_2d.m, triangle_2d.m, box_2d.m, circ_2d.m
% Reference: http://en.wikipedia.org/wiki/Lanczos_resampling
% Example: [MatX,MatY]=meshgrid([-10:1:10],[-10:1:10]);
%          F=lanczos_2d(MatX,MatY,2,1,[0 0],1);
%          surface(F);
% Reliable: 2
%--------------------------------------------------------------------------

Def.A       = 2;
Def.Stretch = 1;
Def.Pos0    = [0 0];
Def.Norm    = 1;

if (nargin==2),
   A       = Def.A;
   Stretch = Def.Stretch;
   Pos0    = Def.Pos0;
   Norm    = Def.Norm;
elseif (nargin==3),
   Stretch = Def.Stretch;
   Pos0    = Def.Pos0;
   Norm    = Def.Norm;
elseif (nargin==4),
   Pos0    = Def.Pos0;
   Norm    = Def.Norm;
elseif (nargin==5),
   Norm    = Def.Norm;
elseif (nargin==6),
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

MatR = sqrt((X-X0).^2 + (Y-Y0).^2)./Stretch;
F    = A.*sin(pi.*MatR).*sin(pi.*MatR./A)./(pi.^2.*MatR.^2);

% set singular point to 1:
I = find(isnan(F) | isinf(F));
F(I) = 1;

% set to zero outside A:
I = find(MatR>A);
F(I) = 0;


if (isnan(Norm)),
   % do not normalize
else
   F = Norm.*F./sumnd(F);
end
