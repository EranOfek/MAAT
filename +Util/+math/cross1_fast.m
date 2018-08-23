function C=cross1_fast(A,B)
% Fast version of cross product of two 3-elements vectors
% Package: Util.math
% Description: cross product of two 3-elements vectors. This is a fast
%              version of the cross.m function. This function will work
%              only between two vectors.
% Input  : - First 3-element vector.
%          - Second 3-element vector.
% Output : - Vector of cross product between the two input vectors.
%            The output vector is always a raw vector.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: cross_fast.m
% Example: C=Util.math.cross1_fast([1 2 3],[1 2 1])
% Reliable: 2
%--------------------------------------------------------------------------

C = [A(2).*B(3) - A(3).*B(2), A(3).*B(1) - A(1).*B(3), A(1).*B(2) - A(2).*B(1)];
%C = [A(2,:).*B(3,:) - A(3,:).*B(2,:), A(3,:).*B(1,:) - A(1,:).*B(3,:), A(1,:).*B(2,:) - A(2,:).*B(1,:)];
