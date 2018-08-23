function [Y,DelY]=polyconf_cov(Par,X,Cov)
%------------------------------------------------------------------------------
% polyconf_cov function                                                 FitFun
% Description: Estimate the 1-sigma confidence interval of a polynomial
%              given the polynomial best fit parameters and its covariance
%              matrix. REQUIRE FURTHER TESTING.
% Input  : - Vector of polynomial coef. [Pn, Pn-1, ... P0], as returned
%            by fitgenpoly.m.
%          - X position in which to evaluate the polynomial and its error.
%          - The covariance matrix of the best fit polynomial.
% Output : - Value of the polynomial at X.
%          - 1 sigma error of the polynomial at X.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Jan 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: P=fitgenpoly([1:1:100].',[1:1:100].'+randn(100,1),1,1);
%          [Y,DelY]=polyconf_cov(P.Par,[1 50 100 200]',P.Cov);          
% Reliable: 2
%------------------------------------------------------------------------------

Deg   = length(Par)-1;
SizeX = size(X);
Xv    = X(:);

% for each value in Xv: Xv.^[Deg:-1:0];
Mat  = bsxfun(@power,Xv,[Deg:-1:0]);
N    = length(Xv);
DelY = zeros(N,1);
Y    = zeros(N,1);
for I=1:1:N,
   DelY(I) = sqrt(Mat(I,:)*Cov*Mat(I,:)');
   Y(I)    = polyval(Par,Xv(I));
end

DelY = reshape(DelY,SizeX);
Y    = reshape(Y,SizeX);

