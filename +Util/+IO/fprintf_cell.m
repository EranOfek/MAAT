function String=fprintf_cell(FID,Format,Cell)
% An fprintf command for a cell vector.
% Package: Util.IO
% Description: An fprintf command for a cell vector.
% Input  : - File identifier or string containing file name to open and
%            close at the end.
%            If empty matrix (i.e., []), then return a string containing
%            the cell contents, but donot write a file.
%          - Print format (e.g., '%s\n').
%          - Cell vector to print.
% Output : - A string containing the cell contents.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Util.IO.fprintf_cell('File.dat','%s\n',{'first line','second line'});
% Reliable: 2
%------------------------------------------------------------------------------
if (isempty(FID)==1)
   OnlyString = 1;
else
   OnlyString = 0;
   if (ischar(FID)==1)
      FileName = FID;
      FID = fopen(FileName,'w');
   else
      FileName = [];
   end
end

N = length(Cell);

String = '';
for I=1:1:N
   switch OnlyString
    case 0
       fprintf(FID,Format,Cell{I});
    otherwise
       % do nothing
   end
   String = [String, sprintf(Format,Cell{I})];
end

switch OnlyString
 case 0
    if (isempty(FileName)==0)
       fclose(FID);  
    end
 otherwise
    % do nothing
end
