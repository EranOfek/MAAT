function [Cell,Header]=parse_html_table(Table,TableInd,Convert,RemoveA,Output,MaxRetCounter)
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
% [Cell,Header]=www.parse_html_table('http://ogle.astrouw.edu.pl/ogle3/ews/2008/ews.html',1,'y','y','cm');
% Reliable: 2
%------------------------------------------------------------------------------
Header = {};

DefV.
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (ischar(Table))
    % assume Table is URL
    URL = Table;
    Table = webread(URL);
end


http://www.johnstonsarchive.net/astro/astmoontable.html


Ista = strfind(lower(Table),'<table');
Iend = strfind(lower(Table),'</table');


Tmp = regexp(Table,'<td>|<td nowrap>|<td colspan=','split');
Tmp = Tmp(2:end);
N   = numel(Tmp);
C   = cell(N,1);
K = 0;
for I=1:1:N
    if (~isempty(strfind(Tmp{I},'<td colspan=')))
        ColSpan = str2double(Tmp{I}(1))
        for J=1:1:ColSpan
            K = K + 1;
            C{K} = NaN;
        end
    else
        
        if (~isempty(strfind(Tmp{I},'<a')))
            % skip <a
            A = regexp(Tmp{I},'>','split');
            A = regexp(A{2},'<','split');
            A = A{1};
        else
            Tmp1 = regexp(Tmp{I},'<','split');
            A = Tmp1{1};
        end
        K = K + 1;
        C{K} = A;
    end
end

