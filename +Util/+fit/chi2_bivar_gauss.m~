function [Chi2,F]=chi2_bivar_gauss(X,Y,Par,Par2,Data,Error);
%------------------------------------------------------------------------------
% chi2_bivar_gauss function                                             FitFun
% Description: Calculate the \chi^2 of a bivariate Gaussian with a data
%              and error matrices:
%              Chi2 = sumnd(((bivar_gauss(X,Y,Pars)-Data)/Error)^2))
% Input  : - X
%          - Y
%          - Vector of parameters:
%            [X0, Y0, SigmaX, SigmaY, Rho, Norm, Const],
%            or alternatively a structure (Par) containing some of
%            the parameters: 
%            .X0, .Y0, .SigmaX, .SigmaY, .Rho .Norm, .Const
%            Default for Norm is 1.
%            Default for Const is 0.
%          - Optioanl structure (Par2) containing the rest of the
%            Gaussian parameters:
%            (.X0, .Y0, .SigmaX, .SigmaY, .Rho, .Norm, .Const). 
%            Note that Par2 overide Par.
%          - Data matrix of values to fit at bivar_gauss(X,Y,Pars...).
%          - Errors corresponding to data matrix.
% Output : - Chi2
%          - Value of bivariate Gaussian at [X,Y].
% See also: bivar_gauss.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       May 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [MatX,MatY] = meshgrid([1:1:10],[1:1:10]);
%          Par.X0 = 5; Par.Y0 = 5;
%          Par2.SigmaX = 2; Par2.SigmaY = 2; Par2.Rho = 0;
%          [Chi2,F]=chi2_bivar_gauss(MatX,MatY,Par,Par2,rand(10,10),rand(10,10));
% Reliable: 2
%------------------------------------------------------------------------------

F    = bivar_gauss(X,Y,Par,Par2);
Chi2 = sumnd(((F - Data)./Error).^2);

