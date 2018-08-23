function Bn=bessel_icoef(N,P)
% Calculate the Bessel interpolation coefficiant.
% Package: Util.interp
% Description: Calculate the Bessel interpolation coefficiant.
% Input  : - Order.
%          - Interpolation factor (P=(X-X0)./H).
% Output : - Bessel interpolation coeff.
% Notes  : Used by interp_diff.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------

%B{1} = inline('P - 0.5','P');

% B{1} = @(P)P;
% B{2} = @(P)0.25.*P.*(P - 1);
% B{3} = @(P)(P - 0.5).*P.*(P - 1)./6;
% B{4} = @(P)0.5.*(P + 1).*P.*(P - 1).*(P - 2)./24;
% B{5} = @(P)(P - 0.5).*(P + 1).*P.*(P - 1).*(P - 2)./120;
% B{6} = @(P)0.5.*(P + 2).*(P + 1).*P.*(P - 1).*(P - 2).*(P - 3)./720;


Bn = zeros(length(N),1);
for I=1:1:length(N)
%   Bn(I) = B{N(I)}(P(I));
   switch N(I)
       case 1
           Bn(I) = P(I);
       case 2
           Bn(I) = 0.25.*P(I).*(P(I) - 1);
       case 3
           Bn(I) = (P(I) - 0.5).*P(I).*(P(I) - 1)./6;
       case 4
           Bn(I) = 0.5.*(P(I) + 1).*P(I).*(P(I) - 1).*(P(I) - 2)./24;
       case 5
           Bn(I) = (P(I) - 0.5).*(P(I) + 1).*P(I).*(P(I) - 1).*(P(I) - 2)./120;
       case 6
           Bn(I) = 0.5.*(P(I) + 2).*(P(I) + 1).*P(I).*(P(I) - 1).*(P(I) - 2).*(P(I) - 3)./720;
       otherwise
           error('Illegal order');
   end

end

% Ieven = find(0.5.*N==floor(0.5.*N));
% Iodd  = find(0.5.*N~=floor(0.5.*N));
% 
% Bn    = zeros(size(N)); 
% % if N is odd
% Bn(Iodd) = (P(Iodd) - 0.5)./N(Iodd) .* factorial(P(Iodd) + 0.5.*N(Iodd) - 1.5) ./ (factorial(P(Iodd) + 0.5.*N(Iodd) - 1.5 - N(Iodd) + 1) .* factorial(N(Iodd) - 1));
% 
% % if N is even
% Bn(Ieven) = 0.5.*factorial(P(Ieven) + 0.5.*N(Ieven) - 1) ./(factorial(P(Ieven) + 0.5.*N(Ieven) - 1 - N(Ieven)) .* factorial(N(Ieven)));

