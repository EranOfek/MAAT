function [Num]=str2num_nan(Str)
%--------------------------------------------------------------------------
% str2num_nan function                                             General
% Description: Convert string to number, and return NaN if not a number
%              or empty.
% Input  : - String or cell array of strings.
% Output : - Number or matrix of numbers (if input is a cell array).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jun 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Bug    : If input contains a number followed by a string (e.g., '12a'),
%          then the string part will be ignored.
% Example: str2num_nan({[],1,'11.2';'a','b','15.1'})
%          str2num_nan('166.1')
%          str2num_nan('A')
% Reliable: 2
%--------------------------------------------------------------------------

if (~iscell(Str))
    Str = {Str};
end

Nstr = numel(Str);
Num  = zeros(size(Str));
for Istr=1:1:Nstr
   if (isempty(Str{Istr}))
      Num(Istr) = NaN;
   else
      if (isnumeric(Str{Istr})==1)
         Num(Istr) = Str{Istr};
      else
          
         TmpNum = sscanf(Str{Istr},'%f');
         if (ischar(TmpNum) || isempty(TmpNum))
            Num(Istr)  = NaN;
         else
            Num(Istr)  = TmpNum;
         end
      end
   end
end