function Flag=isnan_cell(Cell)
% Check if cell elemt is NaN.
% Package: Util.cell
% Description: Given a cell array, return a matrix of flags indicating
%              if each one of the cells is nan.
% Input  : - Cell array.
% Output : - Matrix of flags indicating if each one of the cells is
%            NaN. 1 if NaN, and 0 if not NaN.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Jun 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%-------------------------------------------------------------------------

Inotchar = find(cellfun(@ischar,Cell)==0);
FlagSub = cellfun(@isnan,Cell(Inotchar));

Flag = false(size(Cell));
Flag(Inotchar) = FlagSub;

% old code
%Size = size(Cell);
%Flag = zeros(Size);
%N    = prod(Size);
%
%Flag = zeros(N,1);
%for I=1:1:N,
%   if (isnumeric(Cell{I})),
%      Flag(I) = isnan(Cell{I}(1));
%   else
%      Flag(I) = 0;
%   end
%end
%
%Flag = reshape(Flag,Size);
