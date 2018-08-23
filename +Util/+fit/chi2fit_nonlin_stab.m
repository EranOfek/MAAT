function [Chi2,Resid]=chi2fit_nonlin_stab(X,Fun,Data,DataErr,varargin);
%-----------------------------------------------------------------------------
% chi2fit_nonlin_stab function                                         FitFun
% Description: Stab function for chi2fit_nonlin.m
% Input  : - Vector of free parameters.
%          - Model function name or function handle.
%            Ymodel = Fun(X,varargin);
%          - Vector of observations.
%          - Vector of errors in observations.
%            If scalar, then use the same error for all observations.
%            If empty matrix then use default. Default is 1.
%          * Arbitrary number of additional arguments to pass to Fun.
% Output : - The \chi^2
%          - Vector of residuals.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                      Jan 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%-----------------------------------------------------------------------------

Ycalc = feval(Fun,Data,varargin{:});
Resid = Data - Ycalc;
Chi2  = sum((Resid./DataErr).^2);
