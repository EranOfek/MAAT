function [IsEqual,IsEqualToFirst]=cell_equal(Cell,DelSpace)
% Compare the content of a cell array to its first cell.
% Package: Util.cell
% Description: Given a cell array of strings or numbers, compare the value
%              of all the cells to the first cell and check if they are
%              equal.
% Input  : - A cell array
%          - Delete spaces from strings prior to comparison {true|false}.
%            Default is true.
% Output : - A flag indicating of all the cell elements are identical.
%          - A vector indicating if each cell element is equal to the
%            elelemnt in the first cell.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: cell_equal({'A','A','A'}); cell_equal({1,'2','E'});
% Reliable: 2
%--------------------------------------------------------------------------
import Util.string.*

if (nargin<2)
    DelSpace = true;
end


Ncell = numel(Cell);
Cell1 = Cell{1};
IsEqualToFirst    = false(Ncell,1);
IsEqualToFirst(1) = true;
for Icell=2:1:Ncell
    if (ischar(Cell1) && ischar(Cell{Icell}))
        if (DelSpace)
            IsEqualToFirst(Icell) = strcmp(spacedel(Cell1),spacedel(Cell{Icell}));
        else
            IsEqualToFirst(Icell) = strcmp(Cell1,Cell{Icell});
        end
    elseif (isnumeric(Cell1) && isnumeric(Cell{Icell}))
        if (isnan(Cell1) && isnan(Cell{Icell}))
            IsEqualToFirst(Icell) = true;
        else
            IsEqualToFirst(Icell) = Cell1==Cell{Icell};
        end
    else
        % do nothing - already set to false
    end
end
IsEqual = all(IsEqualToFirst);

    
        
            
