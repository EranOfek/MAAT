function NewCell=ind_cell(Cell,Ind)
% Select indices of vectors in cell array of vectors.
% Package: Util.cell
% Description: Given a cell vector in which each element contains a
%              vector of the same length and a vecor of indices, return
%              a new cell array of the same size in which each element
%              contains a vecor of only the elements which indices are
%              specified (see example).
% Input  : - A cell vecor.
%          - A vector of indices.
% Output : - A new cell vector, in which each element
%              contains a vecor of only the elements which indices are
%              specified (see example).
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Feb 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: A{1} = rand(100,1); A{2} = rand(100,1);
%          I = find(A{1}<0.5);
%          B = ind_cell(A,I);  % returns B{1}=A{1}(I); B{2}=A{2}(I);
% Reliable: 1
%--------------------------------------------------------------------------

Ncell = length(Cell);

NewCell = cell(Ncell,1);
for Icell=1:1:Ncell
    NewCell{Icell} = Cell{Icell}(Ind);
end



