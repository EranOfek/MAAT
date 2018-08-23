function IndVal=bin_sear(X,Val)
% Binary search for a value in a sorted vector.
% Package: Util.find
% Description: Binary search for a value in a sorted vector.
%              If the value does not exist, return the closest index.
% Input  : - sorted vector (ascending).
%          - value to search.
% Output : - Index of closest value.
% See also: find_bin.c
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Sep 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: IndVal=bin_sear([1:1:12]',5.4);
% Reliable: 1
%--------------------------------------------------------------------------
N      = length(X);

if (N==1)
   IndVal = 1;
else
   Ind1   = 1;
   Ind2   = N;
   IndM   = floor(0.5.*N);
   Y1     = X(Ind1);
   Y2     = X(Ind2);
   Ym     = X(IndM);

   Found  = 0;
   while (Found==0)
      %[[Ind1, IndM, Ind2]-499000, [Y1,   Ym,   Y2]]
      if (Val>Ym)
         Ind1 = IndM;
         %Ind2 = Ind2;
         Y1   = X(Ind1);
         %Y2   = X(Ind2);

         if ((Ind2-Ind1)>=2)
            IndM = floor(floor(0.5.*(Ind2+Ind1)));
         else
            Found = 1;
            if (abs(Val-Y1)<abs(Val-Y2))
               IndVal = Ind1;
            else
               IndVal = Ind2;
            end
         end

         Ym   = X(IndM);

      elseif (Val<Ym)
         Ind2 = IndM;
         %Ind1 = Ind1;
         %Y1   = X(Ind1);
         Y2   = X(Ind2);

         if ((Ind2-Ind1)>=2)
            IndM = floor(floor(0.5.*(Ind2+Ind1)));
         else
            Found = 1;
            if (abs(Val-Y1)<abs(Val-Y2))
               IndVal = Ind1;
            else
               IndVal = Ind2;
            end
         end

         Ym   = X(IndM);

      else
         Found  = 1;
         IndVal = IndM;
      end
   end
end


