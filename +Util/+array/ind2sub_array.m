function [ sub ] = ind2sub_array (siz, ndx)
%ind2sub_array Multiple subscripts from linear index.
%   IND2SUB is used to determine the equivalent subscript values
%   corresponding to a given single index into an array.
%
%   ind2sub_array is identical to matlab function IND2SUB except for the
%   fact its output is an array contains the subindices rather than
%   independent variables as the original function does.
%
%   Class support for input IND:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also IND2SUB.
 

siz = double(siz);
lensiz = length(siz);

sub = zeros(length(ndx),lensiz);
    
if lensiz > 2
    k = cumprod(siz);
    for i = lensiz:-1:3,
        vi = rem(ndx-1, k(i-1)) + 1;
        vj = (ndx - vi)/k(i-1) + 1;
        sub(:,i) = double(vj);
        ndx = vi;
    end
end

if lensiz >= 2
    vi = rem(ndx-1, siz(1)) + 1;
    sub(:,2) = double((ndx - vi)/siz(1) + 1);
    sub(:,1) = double(vi);
else 
    sub(:,1) = double(ndx);
end

end

