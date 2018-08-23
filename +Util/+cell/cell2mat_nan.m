function Mat=cell2mat_nan(Cell)
% Convert numeric cell to matrix. Replace empty cells with NaNs.
% Package: Util.cell
% Description: Convert a numeric only cell array to matrix. Replace
%              empty cells by NaN.
% Input  : - Cell array.
% Output : - Matrix.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jan 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%------------------------------------------------------------------------

Size = size(Cell);
Mat  = zeros(Size);

Ie = find(Util.cell.isempty_cell(Cell)==1);
for I=1:1:length(Ie)
   Cell{Ie(I)} = NaN;
end

Mat = cell2mat(Cell);
