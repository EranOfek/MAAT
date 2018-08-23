function [I,Ifound]=findmany(Mat,FindVal)
% Find all values in a vector in another vector or matrix.
% Package: Util.array
% Description: Find all values in a vector in another vector or matrix.
% Input  : - Matrix.
%          - Vector of values. These values will be searched in the matrix.
% Output : - Indices of all the values found in the matrix.
%          - Cell array in which the i-th element is a vector of indices in
%            the matrix that are equal to one value of the i-th element
%            in the vector.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Nov 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Mat = [1;2;3;3;4;4;5];
%          [I,Ifound]=findmany(Mat,[1;4]); 
% Reliable: 1
%--------------------------------------------------------------------------

Nval = length(FindVal);
I    = [];
Ifound = cell(Nval,1);
for Ival=1:1:Nval,
   Ifound{Ival} = find(Mat==FindVal(Ival));
   I      = [I; Ifound{Ival}];
end
