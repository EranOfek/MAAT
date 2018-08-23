function [SortedCell,SortedInd]=sort_numeric_cell(Cell,varargin
% sort each row or columns in a cell array of numbers.
% Package: Util.cell
% Description: sort each row or columns in a cell array of numbers.
%              see also: sortrows_numeric_cell.m
% Input  : - Cell array containing numbers.
%          - Dimension along which to sort (default is 1).
%          - Mode ('ascend' | 'descend'). Default is 'ascend'.
% Output : - A sorted cell array.
%          - The indices such that SortedCell=Cell(SortedInd).
% See also: sortrows_numeric_cell.m
% Tested : Matlab 2012a
%     By : Eran O. Ofek                    Oct 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/ 
% Example: [SortedCell,SortedInd]=Util.cell.sort_numeric_cell(Cell)
% Reliable: 2
%--------------------------------------------------------------------------


Mat = cell2mat(Cell);
[SortedMat,SortedInd] = sort(Mat,varargin{:});
SortedCell = num2cell(SortedMat);



