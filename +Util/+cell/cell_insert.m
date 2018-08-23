function NewCell=cell_insert(Cell,Insert,Ind)
% Insert elements within a cell array.
% Package: Util.cell
% Description: Insert elements within a cell array.
% Input  : - A 1-D cell array.
%          - An element or a cell array of multiple elements to insert.
%            Note that if multiple elements that they will be inserted
%            together at the same position.
%          - Index in which to insert. Use 0 to insert in the begining.
%            If multiple indices, then the whole insert string will be
%            inserted in each requested position.
% Output : - A 1-D cell array with the inserted elements.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: NewCell=cell_insert({'a',b','c','d',e'},{'b1','b2'},2)
%          NewCell=cell_insert({'a',b','c','d',e'},{'b1'},[2 2 2])
% Reliable: 2
%--------------------------------------------------------------------------


if (~iscell(Insert)),
    Insert = {Insert};
end

Nind = numel(Ind);
NewCell = Cell;
for Iind=1:1:Nind,
    N = numel(Cell);
    if (Ind(Iind)==0),
        NewCell = [Insert, NewCell];
    elseif (Ind(Iind)>=N),
        NewCell = [NewCell, Insert];
    else
        NewCell = [NewCell(1:Ind(Iind)), Insert, NewCell(Ind(Iind)+1:end)];
    end
end
