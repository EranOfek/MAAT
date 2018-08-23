function V=mat2vec(M)
% Convert matrix to vector. Use (:) instead.
% Package: Util.array
% Description: Convert matrix to vector.
% Input  : - Matrix.
% Output : - Column Vector.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                       May 2000   
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html 
% Reliable: 1
%-----------------------------------------------------------------------------
V = M(:);

%[SizeI,SizeJ] = size(M);
%V = reshape(M,[SizeI.*SizeJ, 1]);

%V = zeros(SizeI.*SizeJ,1);
%
%for K=1:1:SizeJ,
%   V((K-1).*SizeI+1:K.*SizeI) = M(:,K);
%end

