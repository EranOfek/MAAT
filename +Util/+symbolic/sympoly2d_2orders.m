function [OutOrder,OutCoef,ChPoly]=sympoly2d_2orders(Poly2D,X,Y)
% Convert a 2D symbolic polynomials into vectors of orders and coef.
% Package: Util.symbolic
% Description: Convert a 2D symbolic polynomials into vectors of orders
%              and coef.
% Input  : - A symbolic 2D polynomial.
%          - Fisrt symbolic variable. Default is X.
%          - Second symbolic variable. Default is Y'.
% Output : - A two rows matrix with the orders of the polynomilas.
%            The first line indicate the orders of the X polynomials.
%            The second line indicate the orders of the Y polynomials,
%            corresponding to the X polynomilas.
%          - Vector of coeficiants multiplying the polynomial orders.
%          - Vector of symbolic polynomial childrens.
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: syms X Y;
%          Poly = 0.1*X^2*Y + 2*X*Y + 3*Y;
%          [OutOrder,OutCoef,ChPoly]=Util.symbolic.sympoly2d_2orders(Poly)
% Reliable: 2

if (nargin<2)
    syms X Y;
end


% The component of the polynomial
% ChPoly = children(Poly2D);
ChPoly = children(sum(expand(Poly2D)));
Nch    = numel(ChPoly);  % number of components

OutOrder = zeros(2,Nch);
OutCoef  = zeros(1,Nch);
for Ich=1:1:Nch
    % for each component
    FactorChPoly = factor(ChPoly(Ich));
    %FactorChPoly = collect(ChPoly(Ich));
    
    IX     = has(FactorChPoly,X);
    IY     = has(FactorChPoly,Y);
    Icoef  = ~(IX | IY);
    
    PowerX = sum(IX);
    PowerY = sum(IY);
    
    %ReNorm = NormXY.^PowerX .* NormXY.^PowerY;
    OutCoef(Ich)    = prod(FactorChPoly(Icoef));
    OutOrder(:,Ich) = [PowerX; PowerY];
    
    
end
