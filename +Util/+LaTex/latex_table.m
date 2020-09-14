function latex_table(OutFileName,CellArray,Header,Format,Title,Caption,ColPos)
% Create a latex table from a data given in a cell array.
% Package: Util.LaTex
% Description: Create a latex table from a data given in a cell array.
% Input  : - Output file name.
%          - Cell array of data (see colvec2cellarray.m).
%          - Header cell-array, contains header for each column.
%          - cell of strings of format for each column.
%          - Table title, default is empty string.
%          - File name containing the caption text, default is no caption.
%          - Column position (e.g., 'clrc') one character per column.
%            'c' for center; 'l' for left; 'r' for right.
%            Default is string of 'c' (in length equal to number of columns.
% Output : The output is written to the output file name.
% Tested : Matlab 5.3
%     By : Eran O. Ofek          October 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: C=colvec2cellarray([1;2;3],['a';'b';'c']);
%          Header=cell(2,2);
%          Header{1,1} = 'Num'; Header{1,2}='Str'; Header{2,2}='[km]';
%          Format=cell(1,2);
%          Format{1} = '$%7.2f%$'; Format{2}='%s';
%          Util.LaTex.latex_table('a.out',C,Header,Format,'Table 1.1','a.try','lr');
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin==4)
   Title        = '';
   Caption      = [];
   ColPos       = [];
elseif (nargin==5)
   Caption      = [];
   ColPos       = [];
elseif (nargin==6)
   ColPos       = [];   
elseif (nargin==7)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(Title)==1)
   Title = '';
end

SizeTab    = size(CellArray);
SizeHeader = size(Header);
if (SizeTab(2)==SizeHeader(2))
   % checked ok
else
   error('Number of columns should be equal to number of headers column');
end

if (isempty(ColPos)==1)
   for I=1:1:SizeTab(2)
      ColPos = [ColPos, ['c']];
   end
end

if (SizeTab(2)==length(ColPos))
   % checked ok
else
   error('Number of columns should be equal to number of column position characters');
end


FidOut = fopen(OutFileName,'w');

fprintf(FidOut,'\\begin{table}[h] \n');
fprintf(FidOut,'\\centering \n');
fprintf(FidOut,'\\begin{tabular}{%s} \n',ColPos);
fprintf(FidOut,'\\underline{%s}\n',Title);
% write table header
for I=1:1:SizeHeader(1)
   for J=1:1:SizeHeader(2)
      Str = Header{I,J};
      if (isempty(Str)==1)
         Str = '';
      end
      if (J==SizeTab(2))
         % last column
         fprintf(FidOut,'%s \\\\ \n',Str);
      else
         fprintf(FidOut,'%s & ',Str);
      end
   end
end

% write table data
for I=1:1:SizeTab(1)
   for J=1:1:SizeTab(2)
      Temp = CellArray{I,J};
      if (iscell(Temp)==1)
         Temp = Temp{1};
      end
      switch Format{J}
       case '%s'
          if (isnumeric(Temp)==1)
             
             Temp = num2str(Temp);
          end
       otherwise
          % do nothing
      end
      Str  = sprintf(Format{J},Temp);
      
      if (J==SizeTab(2))
         % last column
         fprintf(FidOut,'%s \\\\ \n',Str);
      else
         fprintf(FidOut,'%s & ',Str);
      end
   end
end

fprintf(FidOut,'\\end{tabular} \n');
fprintf(FidOut,'\\label{} \n');
fprintf(FidOut,'\\caption{');
if (isempty(Caption)==1)
   fprintf(FidOut,'} \n');
else
   FidCap = fopen(Caption,'r');
   
   while (feof(FidCap)==0)
      Line = fgetl(FidCap);
      fprintf(FidOut,'%s\n',Line);
   end
   fprintf(FidOut,'}\n',Str);
   
   fclose(FidCap);
end
   
fprintf(FidOut,'\\end{table} \n');

fclose(FidOut);
