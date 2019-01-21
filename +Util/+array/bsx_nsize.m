function [S,varargout]=bsx_nsize(varargin)
% bsx_nsize gets multiple of array sizes (i.e. result of the size function),
% verifies they are candidates for bsxfun operation, and the size of the
% result array of bsxfun operation (S).   
% The function also padding the sizes of the input arrays with ones, in
% case one of the size arrays is shorter (i.e. its last dimensions are singleton).
%
% The function uses array sizes rather than the arrays themself in order to
% reduce the data traffic between the functions.
    lengths = cellfun(@length,varargin);    
    sizemat = ones(nargin,max(lengths));
    for i=1:nargin,sizemat(i,1:lengths(i))=varargin{i};end
    S = max(sizemat,[],1);
    if (sum(sizemat(sizemat~=S)>1)~=0)
        error('Matrix dimension must agree.');
    end
    
    nout = min(max(nargout,1)-1,nargin);
    for i = 1:nout, varargout{i}=sizemat(i,:); end
end