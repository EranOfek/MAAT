function Val=legendre_poly(X,N);
%-----------------------------------------------------------------------------
% legendre_poly function                                               FitFun
% Description: Evaluate legendre polynomials
% Input  : - Values at which to evaluate the polynomials.
%          - Order of the polynomials (scalar).
% Output : - Evaluation of the polynomials.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    August 2010
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------------

Val = hypergeom([-N,N+1],1,0.5.*(1-X));
