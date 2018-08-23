function SortedS=sort_struct(Struct,Field,Col)
% Sort all the elements in a structure by one fields in the structure. 
% Package: Util.struct
% Description: Sort in ascending order all the elements in a
%              structure by one fields in the structure. 
% Input  : - Structure.
%          - String containing structure field which the structure
%            elements will be sorted according to.
%          - Column number[s] is the structure field, on which
%            the columns will be sorted.
% Output : - Sorted structure.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                  January 2008
%    URL : http:/wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

FN       = fieldnames(Struct);
N        = length(FN);
FieldVal = Struct.(Field);
[Sorted,SI] = sortrows(FieldVal);

for I=1:1:N
   if (I>1)
      Struct = SortedS;
   end
   FieldVal = getfield(Struct,FN{I});
   if (isstruct(FieldVal)==1)
      error('This version doesnt work on recursive structures');
   elseif (iscell(FieldVal)==1)
      SortedS = setfield(Struct,FN{I},{FieldVal{[SI],:}}.');
   else
      SortedS = setfield(Struct,FN{I},FieldVal(SI,:));
   end
end
