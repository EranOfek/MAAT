function [IsEqual,IsEqualEach]=cell2_equal(Cell1,Cell2,DelSpace)
% Check if two cell array have an identical content
% Package: Util.cell
% Description: Given two cell array containing string or numeric, compare
%              the content of the two cells element by element.
% Input  : - First cell array.
%          - Secodn cell array.
%          - Delete spaces from strings prior to comparison {true|false}.
%            Default is true.
% Output : - A flag indicating if the two cells are identical.
%          - A vector indicating if each element of the cells is identical.
%            If the cells have different size then return NaN.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [IsEqual,IsEqualEach]=Util.cell.cell2_equal({'2',1},{'2',1})
% Reliable: 2
%--------------------------------------------------------------------------
import Util.string.*

if (nargin<3)
    DelSpace = true;
end

Ncell1 = numel(Cell1);
Ncell2 = numel(Cell2);
IsEqual = false;

if (Ncell1~=Ncell2)
    % Cell1 and Cell2 don't have the same number of elements
    % do nothing - IsEqual already false
    IsEqualEach = NaN;
else
    
    IsEqualEach = false(Ncell1,1);
    
    for Ic=1:1:Ncell1
        if (ischar(Cell1{Ic}) && ischar(Cell2{Ic}))
            if (DelSpace)
                IsEqualEach(Ic) = strcmp(spacedel(Cell1{Ic}),spacedel(Cell2{Ic}));
            else
                IsEqualEach(Ic) = strcmp(Cell1{Ic},Cell2{Ic});
            end
        elseif (isnumeric(Cell1{Ic}) && isnumeric(Cell2{Ic}))
            if (isnan(Cell1{Ic}) && isnan(Cell2{Ic}))
                IsEqualEach(Ic) = true;
            else
                IsEqualEach(Ic) = Cell1{Ic}==Cell2{Ic};
            end
        else
            % do nothing - already set to false
        end
    end
    IsEqual = all(IsEqualEach);
end
    
