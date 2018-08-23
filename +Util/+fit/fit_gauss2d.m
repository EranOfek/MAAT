function [Beta,Chi2,E,O]=fit_gauss2d(MatX,MatY,MatVal,MatErr,InitPar,WithBack)
%--------------------------------------------------------------------------
% fit_gauss2d function                                             General
% Description: Non-linear fitting of a 2D elliptical Gaussian with
%               background.
% Input  : - Matrix of X coordinates.
%          - Matrix of Y coordinates.
%          - Matrix of Z values to fit.
%          - Vector of initial guess parameters [], or [X, Y], or
%            [X Y Sigma Norm Back].
%            Default is empty.
%          - A flag indicating if to fit background. Default is true.
% Output : - Best fit Gaussian parameters:
%            [Normalization, X0, Y0, SigmaX, SigmaY, Rho, Background]
%          - \chi^2
%          * The rest of the output parameters returned by fminsearch.m 
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: fun_gauss2d.m
% Example: [Beta,Chi2,E,O]=fit_gauss2d(MatX,MatY,MatVal,1,X,Y);
% Reliable: 2
%--------------------------------------------------------------------------

Def.MatErr   = ones(size(MatVal));
Def.InitPar  = [];
Def.WithBack = true;

if (nargin==3),
    MatErr    = Def.MatErr;
    InitPar   = Def.InitPar;
    WithBack  = Def.WithBack;
elseif (nargin==4),
    InitPar   = Def.InitPar;
    WithBack  = Def.WithBack;
elseif (nargin==5),
    WithBack  = Def.WithBack;
else
    % do nothing
end

% estimate initial parameters
Back       = min(MatVal(:));
Norm       = sum(MatVal(:)-Back);
[~,MaxInd] = maxnd(MatVal);
X          = MatX(1,MaxInd(2));
Y          = MatY(MaxInd(1),1);
Rho        = 0;
SigmaX     = sqrt(sum((MatVal(:)-Back).*(MatX(:)-X).^2)./sum(MatVal(:)-Back));
SigmaY     = sqrt(sum((MatVal(:)-Back).*(MatY(:)-X).^2)./sum(MatVal(:)-Back));
%[Normalization, X0, Y0, SigmaX, SigmaY, Rho, Background].


switch numel(InitPar)
    case 0
        GuessBeta = [Norm, X,          Y,          SigmaX, SigmaY, Rho, Back];
    case 2
        GuessBeta = [Norm, InitPar(1), InitPar(2), SigmaX, SigmaY, Rho, Back];
    case 5
        GuessBeta = [InitPar(4), InitPar(1), InitPar(2), InitPar(3), InitPar(3), Rho, InitPar(5)];
    otherwise
        error('Unknown InitPar vector');
end

%Back = median(MatVal(:));
%GuessBeta = [3.*Back X Y 1 1 0 Back];
% GuessBeta

if (WithBack),
    [Beta,Chi2,E,O]=fminsearch_chi2([MatX(:),MatY(:)],MatVal(:),MatErr(:),{@fun_gauss2d},GuessBeta);
else
    [Beta,Chi2,E,O]=fminsearch_chi2([MatX(:),MatY(:)],MatVal(:),MatErr(:),{@fun_gauss2d},GuessBeta(1:6));
end
