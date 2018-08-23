function IndVal=bin_sear2(X,Target)
%--------------------------------------------------------------------------
% bin_sear function                                                 FitFun
% Description: Binary search for a value in a sorted vector.
%              If the value does not exist, return the closes index.
% Input  : - sorted vector (ascending).
%	   - value to search.
% Output : - index of closest value.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Sep 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: IndVal=bin_sear([1:1:12]',5.4);
% Reliable: 2
%--------------------------------------------------------------------------
High      = length(X);
Low       = 1;

%IndVal = NaN;
Found  = 0;
%Counter = 0;
while (Found==0),
    %Counter = Counter + 1;
       Middle = floor(Low + (High-Low).*0.5);
       
       if (Target<X(Middle)),
           High = Middle - 1;
       elseif (Target>X(Middle)),
           Low = Middle + 1;
       else
           Found  = 1;
       end
       
       if (Low<=High)
           Found = 1;
       end
       %[Counter, Low, High, Middle]
       
end
IndVal = Middle;



