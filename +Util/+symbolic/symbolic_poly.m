function SymPoly=symbolic_poly(Orders,X)
% Build a symbolic polynomial
% Package: Util.symbolic
% Description: Construct a symbolic polynomial with given orders
% Input  : - Vector of orders.
%          - Symbolic variable. Default is X.
% Output : - Vector of symbolic polynomial orders.
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: SymPoly=Util.symbolic.symbolic_poly([2 3 4],X)
% Reliable: 2

if (nargin<2),
    syms X
end

N = numel(Orders);
for I=1:1:N,
    SymPoly(I) = X.^Orders(I);
end


