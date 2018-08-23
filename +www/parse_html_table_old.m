function [Cell,Header]=parse_html_table_old(Table,TableInd,Convert,RemoveA,Output,MaxRetCounter)
% Parse columns from an HTML table into matlab
% Package: www
% Description: Parse columns from an HTML table into matlab.
%              The program can not parse tables with colspan parameter
%              different than 1.
% Input  : - String containing URL name from which to parse the table.
%            Alternatively, this can be a file identifier that contains
%            the HTML table. In this case the user is responsible for
%            opening and closing the file. Or this can be a cell array
%            in which there is one cell that contains the html text.
%          - Number of table in HTML page. This is useful if more than one
%            table is found within HTML page. Default is 1.
%            If two element vector is given then assume the table is nested.
%            In this case the first number indicate the number of the
%            outer table and the second is the number of the nested table.
%          - Convert all data in table to double {'y','n'}, default is 'n'.
%          - Remove anchor links from the table {'y' | 'n'}, default is 'y'.
%          - Output type:
%            'cm' - cell matrix.
%            'cv' - cell vector of cells (default).
%          - Maximum number of retrieves. Default is 1.
%            If >1, then will try again in case of a failure to access
%            the URL. 
% output : - Cell table.
%          - Cell array containing table columns header.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Jun 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [Cell,Header]=www.parse_html_table_old('http://ogle.astrouw.edu.pl/ogle3/ews/2008/ews.html',1,'y','y','cm');
% Reliable: 2
%------------------------------------------------------------------------------
Header = {};

Def.TableInd = 1;
Def.Convert  = 'n';
Def.RemoveA  = 'y';
Def.Output   = 'cv';
Def.MaxRetCounter = 1;
if (nargin==1)
   TableInd = Def.TableInd;
   Convert  = Def.Convert;
   RemoveA  = Def.RemoveA;
   Output   = Def.Output;
   MaxRetCounter = Def.MaxRetCounter;
elseif (nargin==2)
   Convert  = Def.Convert;
   RemoveA  = Def.RemoveA;
   Output   = Def.Output;
   MaxRetCounter = Def.MaxRetCounter;
elseif (nargin==3)
   RemoveA  = Def.RemoveA;
   Output   = Def.Output;
   MaxRetCounter = Def.MaxRetCounter;
elseif (nargin==4)
   Output   = Def.Output;
   MaxRetCounter = Def.MaxRetCounter;
elseif (nargin==5)
   MaxRetCounter = Def.MaxRetCounter;
elseif (nargin==6)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isnumeric(Table))
   % Table Is FID
   Str = Util.files.file2str(Table,'str');
elseif (ischar(Table))
   % Table is a URL
   Done = 0;
   RetCounter = 0;
   while (Done==0 && RetCounter<MaxRetCounter)
      RetCounter = RetCounter + 1;
      try
         %Str = urlread(Table);
         Str = webread(Table);
         Done = 1;
      catch
         if (MaxRetCounter>1)
            fprintf('  Failed reading URL - try again in 30s\n');
            pause(30);
         end
      end
   end
elseif (iscell(Table))
   % Table is already string
   Str = Table{1};
else
   error('Unknown Table format');
end

if (length(TableInd)>1)
   TableIndNest = TableInd(2);
   TableInd     = TableInd(1);
else
   TableIndNest = 0;
end
TableInd = TableInd + TableIndNest;

TableStart = strfind(lower(Str),'<table');
TableEnd   = strfind(lower(Str),'</table>')+7;

TableStart = TableStart(TableInd);
TableEnd   = TableEnd(TableInd-TableIndNest);

% select table from Str:
Str = Str(TableStart:TableEnd);

%IIcg1 = strfind(Str,'<colgroup>');
%IIcg2 = strfind(Str,'</colgroup>');
%Str = [Str(1:IIcg1-1), Str(IIcg2+11:end)];


% remove <table>, </table>, </tr>, </td>, </th>
Str = regexprep(Str,{'<table.*?>','</table>','</tr>','</td>','</th>'},'','ignorecase');



% Split by <tr>
SplitTR = regexpi(Str,'<tr','split');
Ntr     = length(SplitTR);
Counter = 0;
for Itr=1:1:Ntr
   if (isempty(strfind(lower(SplitTR{Itr}),'<td')))
      % skip
      if (isempty(strfind(lower(SplitTR{Itr}),'<th')))
         % skip
      else
         % parse header
         SplitTH = regexpi(SplitTR{Itr},'<th.*?>','split');
         Header  = SplitTH(2:end);   % remove first cel that should be empty
      end
   else
      Counter = Counter + 1;
      SplitTH = regexpi(SplitTR{Itr},'<td.*?>','split');
      SplitTH = SplitTH(2:end);   % remove first cel that should be empty

      % remove anchors
      switch RemoveA
          case 'y'
             SpliTH = regexprep(SplitTH,{'<a.*?>','</a>'},'','ignorecase');
          otherwise
             % do nothing
      end

      if (Counter==1)
         AllCell = SplitTH;
      else
        AllCell = [AllCell;SplitTH];
      end
   end
end

switch lower(Convert)
 case 'y'
    AllCell = str2double(AllCell);
 otherwise
    % do nothing
end

switch lower(Output)
 case 'cv'
    Ncol = size(AllCell,2);
    for Icol=1:1:Ncol
       Cell{Icol} = AllCell(:,Icol);
   end
 case 'cm'
    Cell = AllCell;
 otherwise
    error('Unknown Output option');
end
