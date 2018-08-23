function [Vector,Dictionary]=cellstr2num_dictionary(VectorCell,Dictionary)
% Given a cell array of strings, repace unique strings by numeric index
% Package: Util.cell
% Description: Concert a cell array of strings to a vector of numeric
%              indices. Each numeric index represent a unique string.
%              Also returns a dictionary that translate the string name to
%              its index.
% Input  : - Cell array of strings.
%          - Optional dictionary to which to append the new dictionary.
% Output : - A vector of indices.
%          - A structure array containing the dictionary.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Vector,Dictionary]=Util.cell.cellstr2num_dictionary(C{3});
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin<2)
    Dictionary = Util.struct.struct_def({'Ind','Name'},0,1);
end

UniqueStr = unique(VectorCell);
Nun = numel(UniqueStr);
Vector = nan(Nun,1);
for Iun=1:1:Nun
    Flag = strcmp(VectorCell,UniqueStr{Iun});
    
    IndExist = find([Dictionary.Ind] == Iun);
    
    Ndic = numel(Dictionary);
    if (isempty(IndExist))
        % not in dictionary - add
        Dictionary(Ndic+1).Ind  = Iun;
        Dictionary(Ndic+1).Name = UniqueStr{Iun};
        Vector(Flag) = Ndic+1;
    else
        % already in dictionary
        Vector(Flag) = IndExist;
    end
end
