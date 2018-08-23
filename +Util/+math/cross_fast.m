function C=cross_fast(A,B)
% Fast cross product of two 3-elements matrices
% Package: Util.math
% Description: cross product of two 3-columns matrices. This is a fast
%              version of the cross.m function. 
% Input  : - First 3-element vector.
%          - Second 3-element vector.
% Output : - Vector of cross product between the two input vectors.
%            The output vector is always a raw vector.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: cross1_fast.m
% Example: C=Util.math.cross_fast(rand(10,3),rand(10,3))
% Reliable: 2
%--------------------------------------------------------------------------

C = [A(:,2).*B(:,3) - A(:,3).*B(:,2), A(:,3).*B(:,1) - A(:,1).*B(:,3), A(:,1).*B(:,2) - A(:,2,:).*B(:,1)];
