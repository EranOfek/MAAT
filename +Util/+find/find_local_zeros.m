function Zeros=find_local_zeros(X,Y,Deg)
%-----------------------------------------------------------------------------
% find_local_zeros function                                            FitFun
% Description: Given table of equally spaced data, use Stirling
%              interpolation formula to find the local zeros of
%              the tabulated data. The program find all the local
%              zeros between X(Deg/2) to X(end-Deg/2), where
%              Deg is the degree of interpolation.
% Input  : - Equally spaced (and asendingly sorted) X.
%          - Y
%          - Order (currently supports only 4th degree), default is 4.
% Output : - List of all zeros:
%            values of [X, 1st derivative, 2nd derivative d^2Y/dX^2]
%            in the points of zeros.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% BUGS   : Fails to find some zeros
% see also: interp_diff.m, find_local_extramum.m
% Reliable: 2
%-----------------------------------------------------------------------------
Check   = 'n';    % check if X is equally spaced and sorted...
if (nargin==2)
   Deg = 4;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end  

if (Deg~=4)
   error('Deg~=4 is not supported');
end


Istart = 0.5.*Deg+2;
Iend   = length(X) - 0.5.*Deg-1;
Xi     = X(Istart:1:Iend);
N      = length(Xi);


[Yi,SPoly,H,X0] = Util.interp.interp_diff(X,Y,Xi,Deg,Check);


Zeros  = [];
If     = 0;
for I=1:1:N-1
   A4 = SPoly(I,1);
   A3 = SPoly(I,2);
   A2 = SPoly(I,3);
   A1 = SPoly(I,4);
   A0 = SPoly(I,5);

   Poly   = [A4 A3 A2 A1 A0];
   PolyD1 = [4.*A4 3.*A3 2.*A2 A1];
   PolyD2 = [12.*A4 6.*A3 2.*A2];    

   Roots  = roots(Poly);
   %Zeros = zeros(0,3);
   for Ir=1:1:length(Roots)
      if (isreal(Roots(Ir))==1)
         if (Roots(Ir)>=-0.5 && Roots(Ir)<=0.5)
            % non-complex root has been found in the -0.5<=p<=0.5 range
            If     = If + 1;
            D1     = polyval(PolyD1,Roots(Ir));
            D2     = polyval(PolyD2,Roots(Ir));
            Zeros(If,:) = [X0(I)+H.*Roots(Ir), D1, D2];             
         end
      end   
   end
end
