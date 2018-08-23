function Flag=isempty_cell(Cell)
% Check if each cell element is empty.
% Package: Util.cell
% Description: Given a cell array, return a matrix of flags indicating
%              if each one of the cells is empty.
% Input  : - Cell array.
% Output : - Matrix of flags indicating if each one of the cells is
%            empty. 1 if empty, and 0 if not empty.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Flag=Util.cell.isempty_cell({1, [], 'a'})
% Reliable: 2
%--------------------------------------------------------------------------

Flag = cellfun(@isempty,Cell);

% Old code:
%Size = size(Cell);
%Flag = zeros(Size);
%N    = prod(Size);
%
%Flag = zeros(N,1);
%for I=1:1:N,
%   Flag(I) = isempty(Cell{I});
%end
%
%Flag = reshape(Flag,Size);
