function Vec=and_mat(Mat,Dim)
% Perform logical and operation between all the columns or rows of a matrix
% Package: Util.array
% Description: Perform logical and operation between all the columns or
%              rows of a matrix.
% Input  : - Matrix of logicals.
%          - Diemension along to do the logical and operation.
%            Default is 1.
% Output : - Vector of logical result.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Mat = [false true; true true; false false]; Vec=and_mat(Mat,2);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1),
    Dim = 1;
end

Vec = prod(single(Mat),Dim)==1;