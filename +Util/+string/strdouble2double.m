function Out=strdouble2double(In);
%-----------------------------------------------------------------------------
% strdouble2double function                                           General
% Description: Convert string, souble or any data type to double.
%              Unlike str2double.m, this script doesn't return NaN if
%              the input is already a double.
% Input  : - String or a number or a cell array of strings or numbers.
% Output : - The input in double format - will convert cell to matrix.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                   October 2010
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------------

if (iscell(In)),
   IndNum = find(cellfun(@isstr,In)==0);
   Out    = str2double(In);
   Out(IndNum) = [In{IndNum}];

else
   if (isnumeric(In)),
      Out = In;
   else
      Out = str2double(In);
   end
end

