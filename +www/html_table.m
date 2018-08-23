function html_table(Output,varargin)
% Create an HTML table
% Package: www
% Description: Given a matlab matrix or cell array create an html
%              page with the matrix in an html table.
% Input  : - Output file name.
%          * Arbitrary pairs of arguments:...keyword,value,...
%            where keyword can be one of the following:
%            'TableCell'  - Table Cell array.
%            'TableMat'   - Table matrix (dont supply not cell and mat).
%            'TableLink'  - Cell array of table links, default is empty.
%            'TableHead'  - Cell vector of table columns header.
%            'PreTable'   - String of pre table text.
%            'PostTable'  - String of pos ttable text.
%            'PageTitle'  - Page title string.
%            'TableBorder'- Table border number, default is 1.
%            'BgColor'    - HTML background color, defaukt is '#ffffff'.
%            'TextColor'  - HTML text color, defaukt is '#000000'.
%            'LinkColor'  - HTML link color, default is '#0000ff'.
%            'VLinkColor' - HTML vlink color, default is '#ff0000'.
%            'BodyStatment'-HTML (extra) body statment, default is ''.
%            'TableStatment'- Table (extra) statement, default is ''.
%            'HeadStatment' - Cell vector of (extra) statment within <th>.
%            'CellStatment' - Cell array of (extra) statment within <td> (cell content).
% Output : null
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

Narg = length(varargin);
% set default values
TableCell   = NaN;
TableMat    = NaN;
TableLink   = {};
TableHead   = NaN;
PreTable    = '';
PostTable   = '';
PageTitle   = '';
TableBorder = 1;
BgColor     = '#ffffff';
TextColor   = '#000000';
LinkColor   = '#0000ff';
VLinkColor  = '#ff0000';
BodyStatment= '';
TableStatment='';
HeadStatment= {};
CellStatment= {};
for Iarg=1:2:Narg-1
   switch varargin{Iarg}
    case 'TableCell'
       TableCell = varargin{Iarg+1};
    case 'TableMat'
       TableMat = varargin{Iarg+1};
    case 'TableLink'
       TableLink = varargin{Iarg+1};
    case 'TableHead'
       TableHead = varargin{Iarg+1};
    case 'PreTable'
       PreTable = varargin{Iarg+1};
    case 'PostTable'
       PostTable = varargin{Iarg+1};
    case 'PageTitle'
       PageTitle = varargin{Iarg+1};
    case 'TableBorder'
       TableBorder = varargin{Iarg+1};
    case 'BgColor'
       BgColor = varargin{Iarg+1};C
    case 'TextColor'
       TextColor = varargin{Iarg+1};
    case 'LinkColor'
       LinkColor = varargin{Iarg+1};
    case 'VLinkColor'
       VLinkColor = varargin{Iarg+1};
    case 'BodyStatment'
       BodyStatment = varargin{Iarg+1};
    case 'TableStatment'
       TableStatment = varargin{Iarg+1};
    case 'HeadStatment'
       HeadStatment = varargin{Iarg+1};
    case 'CellStatment'
       CellStatment = varargin{Iarg+1};
    otherwise
       error('Unkown keyword option');
   end
end


[N1,M1] = size(TableCell);
[N2,M2] = size(TableMat);
N       = max([N1;N2]);
M       = max([M1;M2]);
%--- create cell array of haed statment ---
if (isempty(HeadStatment)==1)
   for J=1:1:M
         HeadStatment{J} = '';
   end
end
%--- create cell array of cell statment ---
if (isempty(CellStatment)==1)
   for I=1:1:N
      for J=1:1:M
         CellStatment{I,J} = '';
      end
   end
end
%--- create cell array of table link ---
if (isempty(TableLink)==1)
   for I=1:1:N
      for J=1:1:M
         TableLink{I,J} = {};
      end
   end
end




FID = fopen(Output,'w');

%--- write HTML openinig ---
fprintf(FID,'<HTML> \n');
fprintf(FID,'<HEAD> \n');
fprintf(FID,'<TITLE> \n');
fprintf(FID,'%s \n',PageTitle);
fprintf(FID,'</TITLE> \n');
fprintf(FID,'</HEAD> \n');
fprintf(FID,'<BODY bgcolor=%s text=%s link%s vlink=%s BodyStatment=%s> \n',BgColor,TextColor,LinkColor,VLinkColor,BodyStatment);

fprintf(FID,'%s \n',PreTable);


fprintf(FID,'<TABLE border=%d %s> \n',TableBorder,TableStatment);


if (iscell(TableHead)==1)
   fprintf(FID,'<tr> \n');
   M= length(TableHead);
   %--- write table header ---
   for J=1:1:M
      fprintf(FID,'   <th %s> \n',HeadStatment{J});
      fprintf(FID,'      %s \n',TableHead{J});   
   end
end

if (isnan(TableMat)==0 || iscell(TableCell)==1)
   %--- write table content ---

   %--- Write TabelCell ---
   if (iscell(TableCell)==1)
      [N,M] = size(TableCell);
      for I=1:1:N
         fprintf(FID,'<tr> \n');
         for J=1:1:M
            fprintf(FID,'   <td %s> \n',CellStatment{I,J});
            CellContent = TableCell{I,J};

            %--- write link header ---
            if (isempty(TableLink{I,J})==0)
               %--- creat link ---
	       fprintf(FID,'      <a target="_top" href="%s">\n',TableLink{I,J});
            end

            %--- write table cell content ---
            if (ischar(CellContent))
               fprintf(FID,'      %s \n',CellContent);
            else
                if (isinteger(CellContent))
                    fprintf(FID,'      %d \n',CellContent);
                elseif (islogical(CellContent))
                    fprintf(FID,'      %d \n',CellContent);
                else
                    fprintf(FID,'      %f \n',CellContent);
                end
            end

            %--- write link footer ---
 	    if (isempty(TableLink{I,J})==0)
               %--- creat link ---
	       fprintf(FID,'      </a>\n');
            end
         end
      end
   end

   %--- Write TableMat ---
   if (isnan(TableMat)==0)
      [N,M] = size(TableMat);

      for I=1:1:N
         fprintf(FID,'<tr> \n');
         for J=1:1:M

            fprintf(FID,'   <td %s> \n',CellStatment{I,J});

            %--- write link header ---
            if (isempty(TableLink{I,J})==0)
               %--- creat link ---
	       fprintf(FID,'      <a target="_top" href="%s">\n',TableLink{I,J});
            end

            %--- write table cell content ---
            fprintf(FID,'      %f \n',TableMat(I,J));

            %--- write link footer ---
 	    if (isempty(TableLink{I,J})==0)
               %--- creat link ---
	       fprintf(FID,'      </a>\n');
            end
         end
      end
   end
end

fprintf(FID,'</TABLE> \n');
fprintf(FID,'%s \n',PostTable);

%--- write html closing ---
fprintf(FID,'</BODY> \n');
fprintf(FID,'</HTML> \n');


fclose(FID);
