function Y=polysval(P,X);
%------------------------------------------------------------------------------
% polysval function                                                     FitFun
% Description: Evaluate multiple polynomials at multiple points. Similar to
%              polyval.m but allows a matrix of polynomials coefficients in
%              which each rows is a different polynomial coefficients.
%              The first column in P corresponds to the coefficient of the
%              highest power and the last column for the power of zero.
% Input  : - Matrix of polynomial coefficients.
%          - Scalar or vector of values in which to evaluate the
%            polynomials. If scalar then will evaluate the different
%            polynomials at the same point. If vector, then the size
%            of this vector must be identical to the number of rows in P.
%            In this case each polynomial will be evaluated at the
%            corresponding X.
% Output : - Evaluate polynomials.
% Tested : Matlab 7.11
%     By : Eran O. Ofek                     April 2011
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

Y = zeros(max(size(P,1),length(X)),1);

Np = size(P,2);

Y = P(:,1);
for I=2:1:Np,
   Y = X.*Y + P(:,I);
end

