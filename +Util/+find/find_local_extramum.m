function Extram=find_local_extramum(X,Y,Deg)
% Use stirling interpolation to find local extramum in vector.
% Package: Util.find
% Description: Given table of equally spaced data, use Stirling
%              interpolation formula to find the local extramums of
%              the tabulated data. The program find all the local
%              extramums between X(Deg/2) to X(end-Deg/2), where
%              Deg is the degree of interpolation.
% Input  : - Equally spaced (and asendingly sorted) X.
%          - Y
%          - Order (currently supports only 4th degree), default is 4.
% Output : - List of all local extramums:
%            values of [X, Y, 2nd derivative d^2Y/dX^2] in the points of 
%            local extramums.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: interp_diff.m, find_local_zeros.m
% Example: X=(1:1:100).'; Y=(X-10).^2; find_local_extramum(X,Y);
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


Extram = zeros(0,3);
If     = 0;
for I=1:1:N
   A4 = SPoly(I,1);
   A3 = SPoly(I,2);
   A2 = SPoly(I,3);
   A1 = SPoly(I,4);
   A0 = SPoly(I,5);

   PolyD1 = [4.*A4 3.*A3 2.*A2 A1];

   if (~any(isnan(PolyD1)))
       Roots  = roots(PolyD1);
       %Extram = zeros(0,3);
       for Ir=1:1:length(Roots)
          if (isreal(Roots(Ir))==1)
             if (Roots(Ir)>=-0.5 && Roots(Ir)<=0.5)
                % non-complex root has been found in the -0.5<=p<=0.5 range
                If     = If + 1; 
                PolyD2 = [12.*A4 6.*A3 2.*A2];
                D2     = polyval(PolyD2,Roots(Ir));
                D0     = polyval(SPoly(I,:),Roots(Ir)); 
                Extram(If,:) = [X0(I)+H.*Roots(Ir), D0, D2];
             end
          end   
       end
   end
end
