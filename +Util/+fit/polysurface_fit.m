function [ParStruct,Resid,Good]=polysurface_fit(X,Y,Z,ErrZ,Order,varargin)
% Fit a surface using a 2-D polynomials.
% Package: @Util.fit
% Description: Fit a surface using a 2-D polynomials.
%              See additional possibilities in @TranClass.
% Input  : - X axis.
%          - Y axis.
%          - Surface to fit Z(X,Y).
%            NaN values will be ignored in fit.
%          - Error in surface Z-value. If empty matrix then set
%            all errors to one (i.e., equal weight).
%          - Degree of polynomial fit [DegX, DegY, DegXY].
%            Default is [3 3 2].
%            This indicates the highest degree of each term.
%            For example [3 3 2] means:
%            Z = a_00 + a_10.*X + a_20.*X.^2 + a_30.*X.^3 +
%                       a_01.*Y + a_02.*X.^2 + a_03.*X.^3 +
%                       a_11.*X.*Y + a_21.*X.^2.*Y + a_12.*X.*Y.^2 +
%                       A_22.*X.^2.*Y.^2;
%          * Arbitrary number of pairs of input arguments:
%            ...,keyword,value,... where keyword is one of the followings:
%            'SigClip'   - [Min Max] sigma clipping, default is [3.5 3.5],
%                          in units of sigma.
%            'MaxNiter'  - Maximum number of sigma clipping iterations,
%                          default is 1 (no sigma clipping).
%            'FunType'   - Basis function type:
%                          'poly'  - normal polynomials (default).
%                          'legendre' - Legendre polynomials.
%                          'chebyshev'- Chebyshev polynomials.
%            'Method'    - Type of least square fitting:
%                          {'chol' | 'orth'}, default is 'chol'.
%                          See lscov.m for details.
% Output : - Structure containing best fit parameters:
%            .Par    - Best fit parameters.
%            .ParErr - Error in best fit parameters.
%            .Npar   - Number of free pramaeters.
%            .Chi2   - \chi^2
%            .Dof    - Degree of freedoms.
%            .XTerm  - The corresponding degree of the X polynomial
%                      for each fitted parameter.
%            .YTerm  - The corresponding degree of the Y polynomial
%                      for each fitted parameter.
%          - Matrix of residuals from best fit (Data - Model).
%          - Matrix of flags indicating if an element was used in the
%            final least square itearion (1) or not (0).
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Aug 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
SizeZ = size(Z);

Def.ErrZ  = ones(SizeZ);
Def.Order = [3 3 2];
if (nargin==3),
   ErrZ   = Def.ErrZ;
   Order  = Def.Order;
elseif (nargin==4),
   Order  = Def.Order;
else
   % do nothing
end

if (isempty(ErrZ)),
   ErrZ  = Def.ErrZ;
end
if (isempty(Order)),
   Order = Def.Order;
end

if (numel(ErrZ)==1),
   ErrZ = ErrZ.*ones(SizeZ);
end

DefV.SigClip   = [3.5 3.5];
DefV.MaxNiter  = 1;
DefV.FunType   = 'poly';
DefV.Method    = 'chol';
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

[MatX, MatY] = meshgrid(X,Y);

% reshape matrices to columns (by stacking it column by column):
Nel     = numel(Z);
VecZ    = reshape(Z,[Nel,1]);
VecErrZ = reshape(ErrZ,[Nel,1]);
VecX    = reshape(MatX,[Nel,1]);
VecY    = reshape(MatY,[Nel,1]);

% build the design matrix - 
Col.XTerm = zeros;
Col.YTerm = zeros;
Counter = 1;
Col.XTerm(Counter) = 0;
Col.YTerm(Counter) = 0;
H = ones(Nel,1);   % constant term (a_00)

%- X axis (a_*0 terms):
for Ix=1:1:Order(1),
   Counter = Counter + 1;
   Col.XTerm(Counter) = Ix;
   Col.YTerm(Counter) = 0;
   H = [H, polydeg_val(VecX,Ix,InPar.FunType)];
end
%- Y axis (a_0* terms):
for Iy=1:1:Order(2),
   Counter = Counter + 1;
   Col.XTerm(Counter) = 0;
   Col.YTerm(Counter) = Iy;
   H = [H, polydeg_val(VecY,Iy,InPar.FunType)];
end
%- X axis (a_*0 terms):
for Ix=1:1:Order(1),
   for Iy =1:1:Order(2),
      Counter = Counter + 1;
      Col.XTerm(Counter) = Ix;
      Col.YTerm(Counter) = Iy;
      H = [H, polydeg_val(VecX,Ix,InPar.FunType).*polydeg_val(VecY,Iy,InPar.FunType)];
   end
end

% solve LSQ problem:
Niter = 0;
FlagGood = ones(Nel,1);   % select all matrix (first iteration)
while (Niter<InPar.MaxNiter)
   Niter = Niter + 1;
   Igood        = find(FlagGood & isnan(VecZ)==0);
   [Par,ParErr] = lscov(H(Igood,:),VecZ(Igood),VecErrZ(Igood).^2,InPar.Method);
   VecResid     = VecZ - H*Par;
   StdResid     = std(VecResid(Igood));
   FlagGood     = (VecResid< (abs(InPar.SigClip(2)).*StdResid) & VecResid> (-abs(InPar.SigClip(1)).*StdResid));
end

ParStruct.Par    = Par;
ParStruct.ParErr = ParErr;
ParStruct.Npar   = size(H,2);
ParStruct.Chi2   = sum((VecResid./VecErrZ).^2);
ParStruct.Dof    = Nel - ParStruct.Npar;
ParStruct.XTerm  = Col.XTerm;
ParStruct.YTerm  = Col.YTerm;

Resid    = reshape(VecResid,SizeZ);
Good     = reshape(FlagGood,SizeZ);



%---------------------------------------
function Val=polydeg_val(X,Deg,FunType);
%---------------------------------------
Def.FunType = 'poly';
if (nargin==2),
   FunType  = Def.FunType;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end

switch lower(FunType)
 case 'poly'
    Val = X.^Deg;
 case 'legendre'
    % Legendre polynomials.
    Val = legendre_poly(X,Deg);
 case 'chebyshev1'
    % Chebyshev polynomials of the first kind
    Val = chebyshev_poly(X,N,1);
 case 'chebyshev2'
    % Chebyshev polynomials of the second kind
    Val = chebyshev_poly(X,N,2);
 otherise
    error('Unknown FunType option');
end




