function [DataTable,ColCell]=parse_html_table_data(TableStr,varargin)
% Parse html table data section into a table
% Package: www
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DataTable,ColCell]=www.parse_html_table_data(TableStr);
% Reliable: 
%--------------------------------------------------------------------------

DefV.Conv2numeric         = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% read header
ColCell  = regexp(TableStr,'<[tT][hH].*?>(.*?)</[tT][hH]>','tokens');
if (isempty(ColCell))
    ColCell  = regexp(TableStr,'<FIELD.*?>.*?<DESCRIPTION>(.*?)</DESCRIPTION>.*?</FIELD>','tokens');
end
ColCell = table2cell(cell2table(ColCell));

% read data
Rows  = regexp(TableStr,'<[tT][rR].*?>(.*?)</[tT][rR]>','tokens');
Cells = regexp(TableStr,'<[tT][dD].*?>(.*?)</[tT][dD]>','tokens');

Ncell = numel(Cells);
Nrows = numel(Rows);

Ncol  = Ncell./Nrows;
if (Ncol~=floor(Ncol))
    error('problem - non integer number of columns');
end

DataTable = cell2table(reshape(Cells,Ncol,Nrows).');

if (InPar.Conv2numeric)
    Cols = DataTable.Properties.VariableNames;
    Ncol = numel(Cols);
    for Icol=1:1:Ncol
        NumCol = str2double(DataTable.(Cols{Icol}));
        if (all(isnan(NumCol)))
            % not numeric
        else    
            DataTable.(Cols{Icol}) = NumCol;
        end
    end
end


