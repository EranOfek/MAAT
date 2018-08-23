function chebyshev2d_2poly(ChebyX,ChebyY,Coef)
% SHORT DESCRIPTION HERE
% Package: Util
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Util.fun.chebyshev2d_2poly(ChebyX,ChebyY,Coef)
% Reliable: 
%--------------------------------------------------------------------------

MaxOrderX = max(ChebyX);
MaxOrderY = max(ChebyY);
MaxOrder  = max(MaxOrderX, MaxOrderY);

syms x

Ncheby = numel(ChebyX);
ChX = zeros(Ncheby,MaxOrder+1);
ChY = zeros(Ncheby,MaxOrder+1);

for Icheby=1:1:Ncheby
    ChX(Icheby,:) = padarray(sym2poly(chebyshevU(ChebyX(Icheby),x)),[0 MaxOrder-ChebyX(Icheby)],'pre');
    ChY(Icheby,:) = padarray(sym2poly(chebyshevU(ChebyY(Icheby),x)),[0 MaxOrder-ChebyY(Icheby)],'pre');
end
% invert poly order to: [1, x, x.^2,...]
ChX = fliplr(ChX);
ChY = fliplr(ChY);

for Ix=1:1:MaxOrderX+1
    for Iy=1:1:MaxOrderY+1
        % ChebyX order, ChebyY order
        [Ix-1, Iy-1]
        Not0 = abs(sign(ChX(:,Ix))).*abs(sign(ChY(:,Iy)));
        sum...
        find(ChX(:,Ix)~=0
        
        