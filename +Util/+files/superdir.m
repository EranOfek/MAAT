function Dir=superdir(MatchStr)
% A 'dir'-like function that can deal with more types of wild cards
% Package: Util.files
% Description: A version of the matlab 'dir' function that can deal with
%              more sophisticated types of wild cards.
%              For example searching for: 'l*00[19-21].fits'.
% Input  : - String of filenames to match
% Output : - A dir structure (see dir.m)
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    May 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------

Len = length(MatchStr);
Ib1 = strfind(MatchStr,'[');
Ib2 = strfind(MatchStr,']');

Nib1= length(Ib1);
Nib2= length(Ib2);

Dir(1).name    = '';
Dir(1).folder  = '';
Dir(1).date    = '';
Dir(1).bytes   = 0;
Dir(1).isdir   = 0;
Dir(1).datenum = 0;

if (Nib1==0 && Nib2==0)
   FileName = MatchStr;

   Dir1 = dir(FileName);
   N1   = length(Dir1);
   
   if (verLessThan('matlab','R2016b'))
       % use this for old MATALB versions:
       Dir  = Util.struct.struct_def({'name','date','bytes','isdir','datenum'},N1,1);
   else
       Dir  = Util.struct.struct_def({'name','folder','date','bytes','isdir','datenum'},N1,1);
   end
   
   for I1=1:1:N1
	    %Dir(end+1) = deal(Dir1(I1));
        Dir(I1) = Dir1(I1);
   end

elseif (Nib1==1 && Nib2==1)
   I  = 1;
   Im = strfind(MatchStr(Ib1(I):Ib2(I)),'-');
   Start = str2double(MatchStr(Ib1(I)+1:Ib1(I)+Im-2));
   Stop  = str2double(MatchStr(Ib1(I)+Im:Ib2(I)-1));

   for If=Start:1:Stop
      if (Ib1==1)
         FileName1 = '';
      else
         FileName1 = MatchStr(1:Ib1-1);
      end
      if (Ib2==Len)
         FileName2 = '';
      else
         FileName2 = MatchStr(Ib2+1:end);
      end
   
      FormatStr = sprintf('%%s%%0%dd%%s',ceil(log10(Stop)));
      FileName  = sprintf(FormatStr,FileName1,If,FileName2);

      Dir1 = dir(FileName);
      N1   = length(Dir1);
   %Dir  = Util.struct.struct_def({'name','folder','date','bytes','isdir','datenum'},Stop-Start+1,1);
      for I1=1:1:N1
         Dir(end+1) = Dir1(I1);
      end

   end
   Dir = Dir(2:end);
else
   error('Illegal filename string');
end


