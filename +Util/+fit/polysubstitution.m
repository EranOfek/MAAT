function PolyF=polysubstitution(Poly,PolyL);
%--------------------------------------------------------------------------
% polysubstitution function                                         FitFun
% Description: Given a polynomial [a_n*X^n+...+a_1*x+a_0] coefficients
%              [a_n, a_n-1,..., a_1, a_0] and a the coefficients of
%              a linear transformation of the form X=A_1*Z+A_0
%              [A_1, A_0], substitute the second polynomial into the first
%              and find the new coefficients.
%              This function is being used by fitgenpoly.m to change the
%              variables in polynomials.
% Input  : - Polynomial coefficients [a_n, a_n-1,..., a_1, a_0].
%          - Linear transformation coefficients [A_1, A_0].
% Output : - New polynomila coefficients [alpha_n,..., alpha_1, alpha_0].
% Tested : Matlab 2011b
%     By : Eran O. Ofek / Iair Arcavi      Jan 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: PolyF=polysubstitution([2 2 0],[1 0]);
% Reliable: 2
%--------------------------------------------------------------------------

N  = length(Poly) - 1;
N2 = length(PolyL);

if (N2~=2),
    error('Second polynomial must be a first order polynomial');
end

PolyF = zeros(size(Poly));
for I=0:1:N,
    IndI = N-I+1;
    
    if (size(Poly,1)==1),
        K = [I:1:N];
    else
        K = [I:1:N].';
    end
    IndK = N-K+1;
    
    PolyF(IndI) = sum(Poly(IndK).*factorial(K)./( factorial(I).*factorial(K-I) ).*PolyL(1).^I.*PolyL(2).^(K-I));
end
