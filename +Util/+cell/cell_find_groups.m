function [Groups]=cell_find_groups(Cell,DelSpace)
% Find all cell lines with identical values, and return indices of groups.
% Package: Util.cell
% Description: Given a cell array, find all lines with identical values.
%              Each such line defines a group.
% Input  : - Cell array which element contains strings or numbers.
%          - Delete spaces from strings prior to comparison {true|false}.
%            Default is true.
% Output : - Structure of groups. Each structure represent a group of rows
%            in the input cell array which have equal values.
%            The structure contains the following fields:
%            .Conntent  - The row of values defines the group.
%            .ptr       - The indices of the rows that belong to the
%                         group.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Groups,Flag]=cell_find_groups({1 2;1 2;'E' 2;'R', 'V';'E',2})
% Reliable: 2
%--------------------------------------------------------------------------



if (nargin<2)
    DelSpace = true;
end


Nrow = size(Cell,1);

Groups     = Util.struct.struct_def({'Content','ptr'},1,1);

Igroup     = 0;
Flag       = false(Nrow,1);
for Irow1=1:1:Nrow
    if (~Flag(Irow1))
        % not in group yet
        Flag(Irow1) = true;
        Igroup = Igroup + 1;
        Groups(Igroup).Content = Cell(Irow1,:);
        Groups(Igroup).ptr = Irow1;
        
        for Irow2=Irow1+1:1:Nrow

            % check if two lines are equal
            if (~Flag(Irow2))
                IsEqual = Util.cell.cell2_equal(Cell(Irow1,:),Cell(Irow2,:),DelSpace);
                if (IsEqual)
                    % found match - add to group
                    Groups(Igroup).ptr     = [Groups(Igroup).ptr; Irow2];
                    Flag(Irow2)            = true;
                end
            end
        end
    end
end

        

