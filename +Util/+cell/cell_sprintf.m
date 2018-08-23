function NewCell=cell_sprintf(CellFormat,CellPar)
% sprintf on content in cell array
% Pacakfe: Util.cell
% Description: sprintf into strings in a cell array.
% Input  : - A cell array of format stings (e.g., 'My Name is %s').
%          - A cell array of parameters to write into the format strings.
%            A line per format string. If only one line is provided than
%            will use it for all format strings.
% Output : - A cell array of output strings.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: NewCell=cell_sprintf({'%s_%s','%s_%s'},{'A','a';'B','b'})
%          NewCell=cell_sprintf({'%s_%s','%s_%s'},{'A','a'})
% Reliable: 2
%--------------------------------------------------------------------------

N = numel(CellFormat);
NewCell = cell(size(CellFormat));

SizeCP  = size(CellPar,1);

for I=1:1:N,
    CI = min(I,SizeCP);
    NewCell{I} = sprintf(CellFormat{I},CellPar{CI,:});
end