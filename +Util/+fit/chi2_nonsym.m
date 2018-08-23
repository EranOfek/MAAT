function [Chi2,ResidSigma]=chi2_nonsym(Y,Yerr,Ymodel);
%------------------------------------------------------------------------------
% chi2_nonsym function                                                  FitFun
% Description: Given measurments with non-symetric Gaussian errors
%              (i.e., upper and lower errors are not equal), calculate
%              the \chi^{2} relative to a given model.
% Input  : - Vector of measurments.
%          - Two column matrix of lower and upper errors on measurments.
%          - Vector of model at the points of measurments.
% Output : - \chi^{2}
%          - Residuals in units of the errors.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                     April 2012
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

if (size(Yerr,2)==1),
   Yerr = repmat(Yerr,1,2);
end
Yerr = abs(Yerr);  % make sure all the errors are positive

Resid  = Y - Ymodel;
%Ineg = find(Resid<0);
Ipos = find(Resid>0);

Err        = Yerr(:,1);   % negative errors
Err(Ipos)  = Yerr(Ipos,2);

ResidSigma = Resid./Err;
Chi2       = sum(ResidSigma.^2);
