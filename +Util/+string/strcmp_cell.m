function Flag=strcmp_cell(CellStr,Template)
%--------------------------------------------------------------------------
% strcmp_cell function                                             General
% Description: Given two cell arrays of strings, check if each one of the
%              strings in the first cell array exist in the second cell
%              array.
% Input  : - A cell array of strings (or a single string).
%            Each one of the strings in this cell array will be compared
%            with all the strings in the second cell array.
%          - A cell array of strings.
% Output : - A vector of flag of the length of first cell array which
%            indicate for each string in the first cell if it exist
%            in the second cell array.
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    May 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:  Flag=strcmp_cell({'A','B','C'},{'A','AA','B','E'})
% Reliable: 2
%--------------------------------------------------------------------------

if (~iscell(CellStr)),
    CellStr = {CellStr};
end
Ns = numel(CellStr);
Flag = false(Ns,1);
for Is=1:1:Ns,
    Flag(Is) = any(strcmp(Template,CellStr{Is}));
end
