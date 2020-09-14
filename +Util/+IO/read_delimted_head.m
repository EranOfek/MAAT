function [Table,Col]=read_delimted_head(File,varargin)
% Read delimited file with header
% Package: Util.IO
% Description: Read a delimited table in which one of the first lines
%              is the header of the delimited table. The program returns
%              a cell array in which each cell is a column in the table
%              and a structure containing the column header names.
%              The default parameters are optimized to deal with some
%              SQL output.
% Input  : - File name.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            The following keywords are available:
%            'HeadLine'    - The number of the line containing the header.
%                            Default is 2.
%            'StartLine'   - The number of the line in which the table begins.
%                            Default is 4.
%            'Delim'       - Delimiter. Default is '|'.
%            'CommentStyle'- Comment style to ignore. Default is '('.
%            'StringHead'  - A cell array which contains strings.
%                            The substrings will be matched with the headers
%                            and if an header contains this substring
%                            than its corresponding column will be treated
%                            as a column of strings. All the other columns
%                            will be read as numbers.
%                            Default is {'fits_'}.
% Output : - Cell array in which is cell is a column in the table.
%          - Structure in which each field is the column name,
%            and the value of the field is the column index.
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Apr 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%--------------------------------------------------------------------------

DefV.HeadLine   = 2;
DefV.StartLine  = 4;
DefV.Delim      = '|';
DefV.StringHead = {'fits_'}; 
DefV.CommentStyle = '(';

InPar  = InArg.populate_keyval(DefV,varargin,mfilename);

FID = fopen(File,'r');
for I=1:1:InPar.HeadLine-1
   fgetl(FID);
end
HeaderLineStr = fgetl(FID);

fclose(FID);
ColNames = strtrim(regexp(HeaderLineStr,sprintf('\\%s',InPar.Delim),'split'));
Ncol     = length(ColNames);
Col      = cell2struct(num2cell(1:1:Ncol),ColNames,2);


Format = cell(Ncol,1);
[Format{Util.cell.isempty_cell(strfind(ColNames,'fits_'))==1}] = deal('%f');
[Format{Util.cell.isempty_cell(strfind(ColNames,'fits_'))==0}] = deal('%s'); 
FormatStr = '';
for Icol=1:1:Ncol
   FormatStr = sprintf('%s %s',FormatStr,Format{Icol});
end

FID = fopen(File,'r');
Table = textscan(FID,FormatStr,'Delimiter',InPar.Delim,'Headerlines',InPar.StartLine-1,'CommentStyle',InPar.CommentStyle);
fclose(FID);
