function NewMat=delete_ind(Mat,Ind,Dim)
% Delete a column/s or row/s from a matrix.
% Package: Util.array
% Description: Delete a column/s or row/s from a matrix.
% Input  : - Matrix.
%          - Indices of rows or columns to remove.
%          - Dimension: 1 - Delete rows; 2 - Delete columns, default is 1.
% Output : - A new matrix with selected column/s or row/s deleted.
% See also: insert_ind.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                  December 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: NewMat=delete_ind(zeros(4,4),[2 3],1)
% Reliable: 1
%---------------------------------------------------------------------------
if (nargin==2)
   Dim = 1;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end


N      = size(Mat,Dim);
NewInd = setdiff((1:1:N),Ind);
switch Dim
 case 1
    NewMat = Mat(NewInd,:);
 case 2
    NewMat = Mat(:,NewInd);
 otherwise
   error(sprintf('%d-dimension is unspported - use only 1/2-d',Dim));
end
