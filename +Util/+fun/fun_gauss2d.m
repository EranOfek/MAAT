function F=fun_gauss2d(Par,XY,Y)
%--------------------------------------------------------------------------
% fun_gauss2d function                                             General
% Description: Calculate bivariate Gaussian in a 2-D grid. This function
%              is appropriate for minimization as the first input
%              argument is the vector of free parameters.
% Input  : - Vector of parameters:
%            [Normalization, X0, Y0, SigmaX, SigmaY, Rho, Background].
%            Background is the only parameter that have default.
%            Background default value is 0.
%          - Two column matrix [X,Y] in which to calculate the 2D Gaussian.
%          - If three arguments are provided then the second is a matrix
%            of X and the third is a matrix of Y.
% Output : - Values of the 2-D Gaussian in the x-y grid.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: bivar_gauss.m, gauss_2d.m, triangle_2d.m
% Example: [MatX,MatY]=meshgrid((0:1:20),(0:1:20));
%          F=fun_gauss2d([1 10 10 2 2 0.5 1],[MatX(:),MatY(:)]);
%          F=fun_gauss2d([1 10 10 2 2 0.5 1],MatX,MatY);
% Reliable: 2
%--------------------------------------------------------------------------



%[Normalization, X0, Y0, SigmaX, SigmaY, Rho, Background]
Norm       = Par(1);
X0         = Par(2);
Y0         = Par(3);
SigmaX     = Par(4);
SigmaY     = Par(5);
Rho        = Par(6);
if (numel(Par)>6),
    Background = Par(7);
else
    Background = 0;
end

if (nargin==2),
    X = XY(:,1);
    Y = XY(:,2);
elseif (nargin==3),
    X = XY;
else
    % do nothing
end

F = Background + Norm./(2.*pi.*SigmaX.*SigmaY.*sqrt(1-Rho.^2)) .* ...
    exp(-1./(2.*(1-Rho.^2)) .* ...
	((X-X0).^2./SigmaX.^2 + ...
	 (Y-Y0).^2./SigmaY.^2 - ...
	 2.*Rho.*(X-X0).*(Y-Y0)./(SigmaX.*SigmaY)));
