function SumBit=sum_bitor(Matrix,Dim)
% A bitor operation on all lines or rows and return a vector ofbit-wise or.
% Package: Util.array
% Description: Given a 2D array of integers, perform a bitor operation
%              on all lines or rows and return a vector ofbit-wise or.
% Input  : - A unsigned integer matrix.
%          - Dimension of bitor summation. Default is 1.
% Output : - A vector of bit wise or between the lines or rows.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Bit=zeros(5,3,'uint16'); Bit(4,2)=1; Bit(4,3)=2; Bit(3,3)=12;
%          SumBit=Util.array.sum_bitor(Bit,2);
% Reliable: 2
%-----------------------------------------------------------------------------


if (nargin==1)
    Dim = 1;
end
if (Dim==2)
    Matrix = Matrix.';
end

Size   = size(Matrix);
SumBit = zeros(1,Size(2),class(Matrix));
for I=1:1:Size(1)
    SumBit = bitor(SumBit,Matrix(I,:));
end


if (Dim==2)
    SumBit = SumBit.';
end

