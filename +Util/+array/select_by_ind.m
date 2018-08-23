function NewMatrix=select_by_ind(Matrix,Ind)
% Select lines by index, return NaN when index is NaN.
% Package: Util.array
% select_by_ind function                                               General
% Description: Select lines from a matrix, given the indices of the lines
%              contains NaNs. Return lines with NaNs from NaNs indices.
% Input  : - Matrix
%          - Vector of indices.
% Output : - Matrix of selected indices. If index is NaN, than the matrix
%            will contain line of NaNs.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                     March 2010
%    URL: http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------------------

Inan = find(isnan(Ind)==1);
IndN = Ind;
IndN(Inan) = 1;
NewMatrix = Matrix(IndN,:);
NewMatrix(Inan,:) = NaN;
