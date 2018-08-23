function NewCell=remove_cell_element(Cell,El2Remove)
% Remove a list of indices from a cell vector.
% Package: Util.cell
% Description: Remove a list of indices from a cell vector.
% Input  : - Cell vector.
%          - Indices of elements to remove from cell vector.
% Output : - New cell vector..
% Tested : Matlab 7.0
%     By : Eran O. Ofek                      June 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cell = {1;2;3}; NewCell=Util.cell.remove_cell_element(Cell,[1 2]);
% Reliable: 1
%-------------------------------------------------------------------------

N2rem   = length(El2Remove);
Ncell   = length(Cell);

NewInd  = setdiff([1:1:Ncell].',El2Remove);
NewCell = Cell(NewInd);

