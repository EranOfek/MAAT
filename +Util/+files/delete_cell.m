function delete_cell(Cell)
% Delete a list of files listed in a cell array.
% Package: Util.files
% Description: Delete a list of files listed in a cell array.
% Input  : - Cell array containing file names to delete.
% Output : *
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: delete({'try.conv','try.param'})
% Reliable: 2
%--------------------------------------------------------------------------

if (ischar(Cell))
    Cell = {Cell};
end

Nc = numel(Cell);
for Ic=1:1:Nc
    delete(Cell{Ic});
end


