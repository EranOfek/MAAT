function Val=chebyshev_poly(X,N,Kind)
%------------------------------------------------------------------------------
% chebyshev_poly function                                              General
% Description: Evaluate Chebyshev polynomials
% Input  : - Values at which to evaluate the polynomials.
%          - Order of the polynomials (scalar).
%          - Kind of Chebyshev polynomial {1 | 2}.
% Output : - Evaluation of the polynomials.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Aug 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Val=chebyshev_poly([6;5],4,1);
% Reliable: 2
%------------------------------------------------------------------------------
switch Kind
 case 1
    % Tn - first kind:
    Val = 0.5.*( (X-sqrt(X.^2-1)).^N + (X+sqrt(X.^2-1)).^N );
 case 2
    % Un - second kind
    Val = ((X + sqrt(X.^2-1)).^(N+1) - (X - sqrt(X.^2-1)).^(N+1))./(2.*sqrt(X.^2-1));
 otherwise
 error('Unknown Kind for Chebyshev polynomials');
end

