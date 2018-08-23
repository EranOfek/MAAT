function [Res,FunStr,FunFun]=fit_2d_polysurface(X,Y,Z,ErrZ,varargin)
%-----------------------------------------------------------------------------
% fit_2d_polysurface function                                          FitFun
% Description: Fit a 2-D polynomial surface to z(x,y) data.
%              e.g., Z= a + b.*X + c.*X.^3 + d.*X.*Y + e.*Y.^2.
% Input  : - Vector of X.
%          - Vector of Y.
%          - Vector of Z measurments in position (X,Y).
%          - Error in Z. If scalar then use equal weitys. Default is 1.
%            If empty use default.
%          * Arbitrary number of pairs of ...,key,val,...
%            The following keywords are available:
%            'Flag' - Flag per each measurment indicating if to use
%                     the measurment in the fit (1) or not.
%                     Default is vecor of ones of the length of X.
%                     If empty use default.
%            'C'    - Add a constant term in the fit (i.e., Z=Const +...)
%                     {'y'|'n'}. Ddefault is 'y'.
%            'X'    - Row vector of the orders of X polynomials in the fit.
%                     For example: [1 3], will add a term Px1.*X + Px2.*X.^3.
%                     Default is [1 2].
%            'Y'    - Row vector of the orders of Y polynomials in the fit.
%                     Default is [1 2].
%            'XY'   - A two column matix of the cross terms orders.
%                     For example: [1 1; 3 4], will add a term:
%                     Pxy1.*X.*Y + Pxy2.*X.^3.*Y.^4.
%                     Default is [1 1].
%            'MaxIter'- Maximum number of sigma clipping iterations.
%                     Default is 0 (i.e., no sigma clipping).
%            'Method'- Sigma clipping method (see clip_resid.m for details).
%                     Default is 'StD'.
%            'Mean' - Method by which to calculate the "mean" of the sample
%                     in the sigma clipping process
%                     (see clip_resid.m for details). Default is 'median'.
%            'Clip' - Two elements vector containing the lower and upper
%                     values for the sigma clipping. This is [Lower, Upper]
%                     number of sigmas (positive) below/above the mean
%                     (see clip_resid.m for details).
%                     See clip_resid.m for defaults.
%            'StdZ' - Add epsilon to sigma clipping std {'y'|'n'},
%                     default is 'y' (see clip_resid.m for details).
% Output : - Structure with the following fields:
%            .PredZ - Z value predicted by best fit at X,Y.
%            .Resid - Vector of all residuals (including those which
%                     were clipped.
%            .Chi2  - \chi^2 of best fit.
%            .Neq   - number of equations.
%            .Npar  - number of free parameters in the fit.
%            .Dof   - number of degrees of freedom.
%            .Igood - Index of values (Z) that were used in the last fit.
%                     (i.e., were not clipped).
%          - String of fitted function, were P is the vector of free
%            parameters.
%          - Function of fitted function.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Aug 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=rand(100,1); Y=rand(100,1);
%          Z=0.3+0.1.*X+0.01.*X.^2+Y+0.001.*X.*Y;
%          [Res,FunStr]=fit_2d_polysurface(X,Y,Z);
% Reliable: 2
%-----------------------------------------------------------------------------

Def.ErrZ = [];
if (nargin==3),
   ErrZ = Def.ErrZ;
end

DefV.Flag     = [];
DefV.C        = 'y';
DefV.X        = [1 2];
DefV.Y        = [1 2];
DefV.XY       = [1 1];
DefV.MaxIter  = 0;
DefV.Method   = 'StD';
DefV.Mean     = 'Median';
DefV.Clip     = [];
DefV.StdZ     = 'y';
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});  


if (isempty(ErrZ)),
   ErrZ = 1;
end
if (numel(ErrZ)==1),
   ErrZ = ErrZ.*ones(size(Z));
end


if (isempty(InPar.Flag)),
    InPar.Flag = true(size(X));
end


Igood   = find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & InPar.Flag);
N       = length(X);

% build the design matrix
% constant term
FunCounter = 0;
switch lower(InPar.C)
 case 'y'
    Hc = ones(N,1);
    FunCounter = 1;
    FunStr = 'P(1)';
 case 'n'
    Hc = ones(N,0);
    FunCounter = 0;
    FunStr = '';
 otherwise
    error('Unknown C option');
end

% X terms
if (isempty(InPar.X)),
    Hx = zeros(N,0);
else
    Hx = bsxfun(@power,X,InPar.X);
end
% Y terms
if (isempty(InPar.Y)),
    Hy = zeros(N,0);
else
    Hy = bsxfun(@power,Y,InPar.Y);
end
% XY terms
if (isempty(InPar.XY)),
    Hxy = zeros(N,0);
else
    Hxy = bsxfun(@power,X,InPar.XY(:,1).').*bsxfun(@power,Y,InPar.XY(:,2).');
end

H = [Hc, Hx, Hy, Hxy];
[Res.Par,Res.ParErr] = lscov(H(Igood,:),Z(Igood), 1./(ErrZ(Igood).^2) );
Res.PredZ = H*Res.Par;
Res.Resid = Z - Res.PredZ;
Res.Chi2  = sum((Res.Resid(Igood)./ErrZ(Igood)).^2);
Res.Neq   = length(Igood);
Res.Npar  = size(H,2);
Res.Dof   = Res.Neq - Res.Npar;
Res.Igood = Igood;

% sigma clipping
for Iiter=1:1:InPar.MaxIter,
   %Igood = clip_resid(Res.Resid,varargin{:}); % this should also take care of NaNs
   [~,~,FlagGood] = clip_resid(Res.Resid,varargin{:}); % this should also take care of NaNs
   Igood = FlagGood & InPar.Flag;
   [Res.Par,Res.ParErr] = lscov(H(Igood,:),Z(Igood), 1./(ErrZ(Igood).^2) );
   Res.PredZ = H*Res.Par;
   Res.Resid = Z - Res.PredZ;
   Res.Chi2  = sum((Res.Resid(Igood)./ErrZ(Igood)).^2);
   Res.Neq   = length(Igood);
   Res.Dof   = Res.Neq - Res.Npar;
   Res.Igood = Igood;
end

if (nargout>1),
   for Ix=1:1:length(InPar.X),
      FunCounter = FunCounter + 1;
      FunStr = sprintf('%s+%s',FunStr,sprintf('P(%d).*X.^%d',FunCounter,InPar.X(Ix)));
   end
   for Iy=1:1:length(InPar.Y),
      FunCounter = FunCounter + 1;
      FunStr = sprintf('%s+%s',FunStr,sprintf('P(%d).*Y.^%d',FunCounter,InPar.Y(Iy)));
   end
   for Ixy=1:1:size(InPar.XY,1),
      FunCounter = FunCounter + 1;
      FunStr = sprintf('%s+%s',FunStr,sprintf('P(%d).*X.^%d.*Y.^%d',FunCounter,InPar.XY(Ixy,1),InPar.XY(Ixy,2)));
   end
   FunFun = eval(sprintf('@(P,X,Y) %s',FunStr));
end

