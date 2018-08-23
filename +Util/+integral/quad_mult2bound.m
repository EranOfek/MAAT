function Y=quad_mult2bound(Fun,A,B,varargin)
% Numerical integration using quad, where the upper bound is a vector
% Package: Util.integral
% Description: Numerical interation using the built in quad.m function
%              where the upper bound of the interation is a vector.
%              To speed the calculation the program sorts the upper bound
%              and perform the integrations between successive bounderies.
% Input  : - Function to integrate (see quad.m for details).
%          - Scalar of lower bound of integration.
%          - Vector of upper bounds of integration.
%          * Additional parameters (see quad.m for details).
% Output : - Integrals.
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    May 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%----------------------------------------------------------------------------

[SB,SI] = sort(B);
N = length(B);

Y = zeros(N,1);
for I=1:1:N
   if (I==1)
      Y(I) = quad(Fun,A,SB(I),varargin{:});
   else
      Y(I) = quad(Fun,SB(I-1),SB(I),varargin{:});
   end
end
Y = cumsum(Y);
Y = Y(SI);
