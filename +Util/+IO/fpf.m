function UsedFormatStr=fpf(FileName,Data,FormatStr,Control)
% Easy to use fprintf, with automatic formatting.
% Package: Util.IO
% Description: Easy to use fprintf, with automatic formatting. This function
%              is similar to fprintf, but (i) open and close the file
%              automaticaly; (ii) in case that format string
%              is not given then the function try to select a format
%              string automaticaly based on the precision of each number.
% Input  : - File name to write the data into, or file identifier
%            (in this case it is the user responsibility to open 
%             and close the file).
%          - Data matrix.
%          - Format string (e.g., '%8.3f %6.3f %5.3f\n').
%            In case format string is not given or empty, then the format
%            string is constructed by estimating a resnoable format.
%          - Control sequence [MaxDexForFloat, 
%                              MaxDexForInteger,
%                              NumberOfBlanks,
%                              Mentisa,
%                              Digit];
%            where:
%            MaxDexForFloat     - calculating the log10 of the ratio
%                                 between the biggest and smallest
%                                 number (=order). If the 'order is
%                                 larger than MaxDexForFloat use 
%                                 exponent number representation, 
%                                 otherwise use float or integer.
%                                 Default is 4.
%            MaxDexForInteger   - If the column contains integers
%                                 and the 'order' is smaller than
%                                 MaxDexForInteger then use integer
%                                 representation, otherwise use
%                                 exponent.
%                                 Default is 6.
%            NumberOfBlanks     - Number of blanks between columns.
%                                 Default is 1.
%            Mentisa            - Default number of digits in mentisa.
%                                 Default is 8.
%            Digit              - Default number of digits after the
%                                 digit.
%                                 Default is 5.
% Output : - The format string used.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Mar 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------
DefMaxDexForFloat   = 4;
DefMaxDexForInteger = 6;
DefNumberOfBlanks   = 1;
DefMentisa          = 8;
DefDigit            = 5;

if (nargin==2)
   FormatStr = [];
   Control   = [DefMaxDexForFloat, DefMaxDexForInteger, DefNumberOfBlanks, DefMentisa, DefDigit];
elseif (nargin==3)
   Control   = [DefMaxDexForFloat, DefMaxDexForInteger, DefNumberOfBlanks, DefMentisa, DefDigit];
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

MaxDexForFloat   = Control(1);
MaxDexForInteger = Control(2);
NumberOfBlanks   = Control(3);
Mentisa          = Control(4);
Digit            = Control(5);

[Ni,Nj] = size(Data);
if (isempty(FormatStr)==1)
   %--- select FormatStr ---
   UsedFormatStr = [];
   for I=1:1:Nj
      %--- check for min/max ---
      Negative = min(Data(:,I))<0;  %1-negative; 0-positive
      Min      = min(abs(Data(:,I)));
      Max      = max(abs(Data(:,I)));
      Dex      = ceil(log10(Max./Min));
      %--- check for intergers ---
      Integer = (sum(floor(Data(:,I))~=Data(:,I))==0); %1-Integer; 0-Float

      if (Integer==1)
         if (Dex<=MaxDexForInteger)
            %--- integer ---
            Format = 'integer';
         else
            %--- exponent ---
            Format = 'exponent';
         end
      else
         if (Dex<=MaxDexForFloat)
            %--- float ---
            Format = 'float';
         else
            %--- exponent ---
            Format = 'exponent';
         end
      end
      if (abs(floor(log10(Min)))>MaxDexForFloat)
         Format = 'exponent';
      end

      %--- column format ---
      switch Format
       case 'integer'
          ColFormat = sprintf('%%%dd',ceil(log10(Max))+Negative);
       case 'exponent'
          ColFormat = sprintf('%%%d.%de',Mentisa,Digit);
       case 'float'
          D0 = Mentisa + Negative + 2;
          D1 = ceil(log10(Max))+Negative;
          D2 = D0 - D1 - 1;
          ColFormat = sprintf('%%%d.%df',D0,D2);

       otherwise
          error('Unknown Format option');
      end

      UsedFormatStr = [UsedFormatStr,blanks(NumberOfBlanks),ColFormat];      
   end
   UsedFormatStr = [UsedFormatStr,'\n'];
else
   UsedFormatStr = FormatStr;
end



if (ischar(FileName)==1)
   FID = fopen(FileName,'w');
else
   FID = FileName;
end
fprintf(FID,UsedFormatStr,Data.');
if (ischar(FileName)==1)
   fclose(FID);
end
