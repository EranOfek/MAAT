function [X,ErrX,Stat]=chi2fit_nonlin(Fun,Data,DataErr,X0,varargin)
%-----------------------------------------------------------------------------
% chi2fit_nonlin function                                              FitFun
% Description: Perform a non-linear \chi^2 fit to dataset.
% Input  : - Model function name or function handle, Ymodel = Fun(X,varargin)
%          - Vector of observations.
%          - Vector of errors in observations.
%            If scalar, then use the same error for all observations.
%            If empty matrix then use default. Default is 1.
%          - Vector of free parameters initial value.
%          * Arbitrary number of additional arguments to pass to Fun.
% Output : - The best fit parameters.
%          - Errors in best fit parameters.
%          - Structure containing fit statistics information,
%            with the following fields:
%            .Chi2   - Chi2 of best fit.
%            .Dof    - Degree of freedom.
%            .Npar   - Number of free parameters.
%            .Nobs   - Number of observations.
%            .Resid  - Vector of residuals.
%            .Hess   - Hessian matrix.
%            .Cov    - Covariance matrix.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                       Jan 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%-----------------------------------------------------------------------------
Def.DataErr = 1;

if (isempty(DataErr)),
   DataErr     = Def.DataErr;
end

if (length(DataErr)==1),
   DataErr = DataErr.*ones(size(Data));
end

X           = fminsearch(@chi2fit_nonlin_stab,X0,[],Fun,Data,DataErr,varargin{:});
X
Stat.Hess   = calc_hessian(@chi2fit_nonlin_stab,X,[],Fun,Data,DataErr,varargin{:});
Stat.Cov    = inv(0.5.*Stat.Hess);
ErrX        = sqrt(diag(Stat.Cov));
[Stat.Chi2,Stat.Resid] = feval(@chi2fit_nonlin_stab,X,Fun,Data,DataErr,varargin{:});
Stat.Nobs   = length(Data);
Stat.Npar   = length(X);
Stat.Dof    = Stat.Nobs - Stat.Npar;



