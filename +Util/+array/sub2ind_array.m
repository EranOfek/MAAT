function ndx = sub2ind_array(siz,v)
%sub2ind_array Linear index from multiple subscripts.
%   sub2ind_array is used to determine the equivalent single index
%   corresponding to a given set of subscript values.
%
%   sub2ind_array is identical to matlab function sub2ind except for the
%   fact its input is an array contains the subindices rather than
%   independent variables as for the original function.
%
%
%   Class support for inputs I,J: 
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also IND2SUB, sub2ind_array.

siz = double(siz);
lensiz = length(siz);
if lensiz < 2
    error(message('MATLAB:sub2ind:InvalidSize'));
end

numOfIndInput = size(v,2);
if lensiz < numOfIndInput
    %Adjust for trailing singleton dimensions
    siz = [siz, ones(1,numOfIndInput-lensiz)];
elseif lensiz > numOfIndInput
    %Adjust for linear indexing on last element
    siz = [siz(1:numOfIndInput-1), prod(siz(numOfIndInput:end))];
end

if any(min(v(:)) < 1) || any(max(v,[],1)>siz)
    %Verify subscripts are within range
    error(message('MATLAB:sub2ind:IndexOutOfRange'));
end

k = [1 cumprod(siz(1:end-1))];
ndx = sum((v-1).*k,2)+1;

% ndx = double(v);
% s = size(v);
% if numOfIndInput >= 2
%     if ~isequal(s,size(v2))
%         %Verify sizes of subscripts
%         error(message('MATLAB:sub2ind:SubscriptVectorSize'));
%     end
%     if any(min(v2(:)) < 1) || any(max(v2(:)) > siz(2))
%         %Verify subscripts are within range
%         error(message('MATLAB:sub2ind:IndexOutOfRange'));
%     end
%     %Compute linear indices
%     ndx = ndx + (double(v2) - 1).*siz(1);
% end 
%     
% if numOfIndInput > 2
%     %Compute linear indices
%     k = cumprod(siz);
%     for i = 3:numOfIndInput
%         v = varargin{i-2};
%         %%Input checking
%         if ~isequal(s,size(v))
%             %Verify sizes of subscripts
%             error(message('MATLAB:sub2ind:SubscriptVectorSize'));
%         end
%         if (any(min(v(:)) < 1)) || (any(max(v(:)) > siz(i)))
%             %Verify subscripts are within range
%             error(message('MATLAB:sub2ind:IndexOutOfRange'));
%         end
%         ndx = ndx + (double(v)-1)*k(i-1);
%     end
% end