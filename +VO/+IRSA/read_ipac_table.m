function [Table]=read_ipac_table(File,FileType,ConvertToJD)
% Read IPAC/IRSA table format
% Package: VO.IRSA
% Description: Read IPAC/IRSA table format from file or string.
% Input  : - A file name containing an IPAC text table.
%          - File type: {'file','str'}. Default is 'file'.
%          - Convert date to JD {true|false}. Default is true.
% Output : - A structure containing the following Fields:
%            'CatCell' - A cell array of the table. Each cell per column.
%            'Col' - Astructure of column indices.
%            'ColCell' - A cell array of column names.
%            'ColType' - A cell array of column types.
%            'ColUnit' - A cell array of column units.
% Reference: http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    May 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Table]=VO.IRSA.read_ipac_table('tt');
% Reliable: 2
%--------------------------------------------------------------------------
import Util.string.*

Def.FileType    = 'file';
Def.ConvertToJD = true;
if (nargin==1)
    FileType    = Def.FileType;
    ConvertToJD = Def.ConvertToJD;
else
    ConvertToJD = Def.ConvertToJD;
end



switch lower(FileType)
    case 'str'
        CellLines = Util.string.strlines2cell(File);
    case 'file'
        % do nothing
        CellLines = Util.files.file2str(File,'cell');
    otherwise
        error('Unknwon FileType option');
end

NotStartTable = true;
Ihead     = 0;
Iheadcol  = 0;
Iline     = 0;
Header    = {};
while NotStartTable
    Iline = Iline + 1;
    Line = CellLines{Iline};
   
    switch Line(1)
        case '\'
            Ihead = Ihead + 1;
            Header{Ihead} = Line(2:end);
        case '|'
            Iheadcol = Iheadcol + 1;
            ColHead{Iheadcol} = Line;
        otherwise
            NotStartTable = false;
    end
    
end

Table.Header = Header;
% construct columns type and description
SpC{1} = regexp(ColHead{1},'\|','split');
SpC{2} = regexp(ColHead{2},'\|','split');
SpC{3} = regexp(ColHead{3},'\|','split');

SpC{1} = spacedel(SpC{1}(2:end-1));
SpC{2} = spacedel(SpC{2}(2:end-1));
SpC{3} = spacedel(SpC{3}(2:end-1));
if (Iheadcol>3)
    SpC{4} = regexp(ColHead{4},'\|','split');
    SpC{4} = spacedel(SpC{4}(2:end-1));
end
Pos    = strfind(ColHead{1},'|');


Ncol = numel(Pos)-1;
%Format = '';
Format = cell(Ncol,3);
for Icol=1:1:Ncol
    switch lower(SpC{2}{Icol})
        case {'r','real','double','float'}
            ColFormat = '%f';
        case {'i','int','long'}
            ColFormat = '%d';
        case {'char','date'}
            ColFormat = '%c';
        otherwise
            error('Unknwon Format type');
    end
    [Format{Icol,1:3}] = deal(Pos(Icol)+1,Pos(Icol+1)-1,ColFormat);
    %Format = sprintf('%s %s',Format,ColFormat);
    Table.Col.(SpC{1}{Icol}) = Icol;
end
%Format = sprintf('%s\\n',Format);

Table.ColCell = SpC{1};
Table.ColType = SpC{2};
Table.ColUnit = SpC{3};

%Table.Cat = Util.string.read_formatted(File,Format,'comment',{'|','\'}); 
Table.Cat = Util.string.read_str_formatted(CellLines,Format,'comment',{'|','\'});

% convert date to JD
if (ConvertToJD)
    IndDate = find(strcmpi(Table.ColType,'date'));
    for Id=1:1:length(IndDate)
        Date = datevec(Table.Cat(:,IndDate(Id)),'yyyy-mm-dd HH:MM:SS');
        Date = Date(:,[3 2 1 4 5 6]);
        JD = celestial.time.julday(Date);
        Table.Cat(:,IndDate(Id)) = JD;
        Table.ColType{IndDate(Id)} = 'double';
        Table.Col.JD = IndDate(Id);
        Table.ColCell{IndDate(Id)} = 'JD';
        Table.ColUnit{IndDate(Id)} = 'day';
    end
end    
