function OutMat=nangetind(Mat,varargin)
% Replace value with NaNs when index is out of bound.
% Package: Util.array
% Description: Get elements from an array by its indices.
%              However, unlike A(I,J) operation in matlab if I or J are
%              out of bound then return NaN.
% Input  : - Column vector, matrix or higher dimension (non-singulation)
%            array.
%          * Arbitrary number of input argument equal to the number
%            of dimensions in the input array.
%            Each contains the indices to retrieve.
% Output : - Retrived array.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Jan 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: A=[1 2 3;4 5 6;7 8 9];
%          OutMat=Util.array.nangetind(A,[-1 2;3 4],[1 2;3 5]);
%          % -1,1 and 4,5 are oput of bounds therefore return: [NaN 5;9 NaN]
%          OutMat=Util.array.nangetind(A,-1,3)
% Reliable: 2
%---------------------------------------------------------------------------


Narg = length(varargin);
Size = size(Mat);

if (Size(1)==1 && Size(2)>1),
   error('If input is vector, then it should be a column vector');
end

AllInd = [];
Ind    = cell(Narg,1);
for Iarg=1:1:Narg,
   Ind{Iarg} = find(varargin{Iarg}>Size(Iarg) | varargin{Iarg}<1);
   AllInd = [AllInd; Ind{Iarg}];
   varargin{Iarg}(Ind{Iarg}) = 1;
end


Ind1   = sub2ind(Size,varargin{:});
OutMat = Mat(Ind1);

OutMat(AllInd) = NaN;


