function F=gauss_2d(X,Y,Sigma,Rho,Pos0,MaxRad,Norm)
%--------------------------------------------------------------------------
% gauss_2d function                                                General
% Description: Calculate bivariate Gaussian in a 2-D grid.
% Input  : - Scalar, vector or matrix of X-coordinates in which to calculate
%            the 2-D Gaussian.
%          - same as the x-ccordinates, but for the y-axis.
%          - Sigma of the Gaussian or [SigmaX, SigmaY] in case sigma
%            is different for each axis.
%            By default SigmaY=SigmaX.
%          - Correlation coef., default is 0.
%            If empty matrix use default.
%          - Center of the Gaussian [X, Y].
%            By default Y=X.
%            Default is [0 0].
%            If empty matrix use default.
%          - Maximum radius of Gaussian behond to set it to zero.
%            Default is Inf.
%            MaxRad is measured from the center of the kernel and not
%            the center of the Gaussian.
%          - A normalization constant of the final matrix.
%            If NaN then donot normalize. Default is 1.
% Output : - Value of the 2-D Gaussian in the x-y grid.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Jul 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: bivar_gauss.m, triangle_2d.m, circ_2d.m, lanczos_2d.m
% Example: [MatX,MatY]=meshgrid([-10:1:10],[-10:1:10]);
%          F=gauss_2d(MatX,MatY,[1],0,[0 0],10);
%          surface(F);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Rho    = 0;
Def.Pos0   = [0 0];
Def.MaxRad = Inf;
Def.Norm   = 1;

if (nargin==3),
   Rho     = Def.Rho;
   Pos0    = Def.Pos0;
   MaxRad  = Def.MaxRad;
   Norm    = Def.Norm;
elseif (nargin==4),
   Pos0    = Def.Pos0;
   MaxRad  = Def.MaxRad;
   Norm    = Def.Norm;
elseif (nargin==5),
   MaxRad  = Def.MaxRad;
   Norm    = Def.Norm;
elseif (nargin==6),
   Norm    = Def.Norm;
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(Pos0)),
   Pos0    = Def.Pos0;
end
if (isempty(MaxRad)),
   MaxRad  = Def.MaxRad;
end
if (isempty(Norm)),
   Norm    = Def.Norm;
end

if (length(Sigma)==1),
   Sigma = [Sigma, Sigma];
end
if (length(Pos0)==1),
   Pos0 = [Pos0, Pos0];
end

SigmaX = Sigma(1);
SigmaY = Sigma(2);

X0 = Pos0(1);
Y0 = Pos0(2);

F = 1./(2.*pi.*SigmaX.*SigmaY.*sqrt(1-Rho.^2)) .* ...
    exp(-1./(2.*(1-Rho.^2)) .* ...
	((X-X0).^2./SigmaX.^2 + ...
	 (Y-Y0).^2./SigmaY.^2 - ...
	 2.*Rho.*(X-X0).*(Y-Y0)./(SigmaX.*SigmaY)));

% set elements outside MaxRad to zero:
if (~isinf(MaxRad)),
   MatR = sqrt(X.^2 + Y.^2);
   I = find(MatR>MaxRad);
   F(I) = 0;
end

if (isnan(Norm)),
   % do not normalize
else
   F = Norm.*F./sumnd(F);
end
