function [ParS,Resid]=fit_gauss2dl(X,Y,Z,ErrZ,FitBack,Method,Elliptical)
%------------------------------------------------------------------------------
% fit_gauss2dl function                                                 FitFun
% Description: Fit a 2-D Gaussian to data of the form f(x,y).
%              f(x,y) = A*exp(-(a*(x-x0)^2+2*b*(x-x0)*(y-y0)+c*(y-y0)^2))
% Input  : - Vector describing the X axis of f(x,y).
%          - Vector describing the Y axis of f(x,y).
%          - Matrix or vector describing f(x,y).
%          - Errors in f(x,y). Default is 1 (equal weight).
%            If empty matrix then use default.
%          - Fit background B + A*exp(...) {'y' | 'n'}, default is 'y'.
%          - Linear least square method:
%            'chol' - Cholesky decomposition (default).
%            'orth' - Orthogonal decomposition.
%          - Fit elliptical gaussian {'y' | 'n'}, default is 'y'.
% Output : - Structure containing the best fit solution.
%            The following fields are available:
%            .Par - Best fit parameters
%                   [A; 1./(2.*(SigmaX.^2 + SigmaY.^2))]
%            .ParErr - Errors in best fit parameters.
%            .Chi2   - \chi^2 of fit.
%            .Dof    - Number of degrees of freedom.
%          - Residuals matrix.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    August 2010
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: % Construct a Gaussian with noise and fit:
%          A=0.01.*randn(100,100);
%          X = [1:1:100]; Y=[1:1:100];  [MatX,MatY]=meshgrid(X,Y);
%          A = A+5.*exp(-((MatX-45).^2 + (MatY-56).^2)./5)
%    [ParS,Resid]=fit_gauss2dl(X,Y,Z,ErrZ,FitBack,Method);
% Reliable: 2
%------------------------------------------------------------------------------

Def.ErrZ       = [];
Def.FitBack    = 'y';
Def.Method     = 'chol';
Def.Elliptical = 'y';
if (nargin==3),
   ErrZ       = Def.ErrZ;
   FitBack    = Def.FitBack;
   Method     = Def.Method;
   Elliptical = Def.Elliptical;
elseif (nargin==4),
   FitBack    = Def.FitBack;
   Method     = Def.Method;
   Elliptical = Def.Elliptical;
elseif (nargin==5),
   Method     = Def.Method;
   Elliptical = Def.Elliptical;
elseif (nargin==6),
   Elliptical = Def.Elliptical;
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end


if (isempty(ErrZ)),
   ErrZ = 1;
end
SizeZ = size(Z);
Nel   = numel(Z);
if (numel(ErrZ)==1),
   ErrZ = ErrZ.*ones(SizeZ);
end

if (sum(SizeZ==1)>0),
   % Z is a vector
   VecX = X;
   VecY = Y;
   VecZ = Z;
   VecErrZ = ErrZ;
else
   % Z is a matrix
   [MatX, MatY] = meshgrid(X,Y);
   VecX = reshape(MatX,[Nel, 1]);
   VecY = reshape(MatY,[Nel, 1]);
   VecZ = reshape(Z,[Nel, 1]);
   VecErrZ = reshape(ErrZ,[Nel, 1]);
end


% Z      = A*exp(-(X.^2+Y.^2)./(2.*(SigmaX.^2 + SigmaY.^2)));
% log(Z) = log(A) -(X.^2+Y.^2)./(2.*(SigmaX.^2 + SigmaY.^2))

% transform to log(f(x,y)):
LogVecErrZ = VecErrZ./VecZ;

switch lower(FitBack)
 case 'n'
    % no background - linear least squares
    LogVecZ = log(VecZ);

    switch lower(Elliptical)
     case 'n'
        % Fit circular Gaussian
        % a*(X-X0)^2 + a*(Y-Y0)^2
        % =
        % X^2 (a)
        % Y^2 (a)
        % X   (-2*a*X0)
        % Y   (-2*a*Y0)
        %     (a*X0^2 + a*Y0^2)
        %---
        %   A+(a*X0^2 + a*Y0^2),    a             ,-2*a*X0,-2*a*Y0
        H = [ones(Nel,1),       -(VecX.^2+VecY.^2), VecX,   VecY];

     case 'y'
        % Fit elliptical Gaussian
        % X^2 (a)
        % Y^2 (c)
        % X   (-2*a*X0 - 2*b*Y0)
        % Y   (-2*b*X0 - 2*c*Y0)
        % X*Y (2*b)
        %     (a*X0^2 + c*Y0^2 + 2*b*X0*Y0)
        %---
        %   
        H = [ones(Nel,1),       -VecX.^2, VecY.^2, VecX,   VecY, VecX.*VecY];

     otherwise
        error('Unknown Elliptical option');
    end
 
    [Par,ParErr] = lscov(H,VecY,1./LogVecErrZ.^2,Method)

    VecResid       = exp(LogVecZ - H*Par);
    ParS.Par       = Par;
%    ParS.Par(1)    = exp(ParS.Par(1));
    ParS.ParErr    = ParErr;
%    ParS.ParErr(1) = ParS.ParErr(1).*exp(ParS.Par(1));
    ParS.Chi2      = sum(VecResid.^2./(VecErrZ.^2));
    Npar           = length(ParS.Par);
    ParS.Dof       = Nel - Npar;
    Resid          = reshape(VecResid,SizeZ);

 case 'y'
    % background - non-linear least squares
    % fminsearch...
    error('FitBack not implemented');
 otherwise
    error('Unknown FitBack option');
end



