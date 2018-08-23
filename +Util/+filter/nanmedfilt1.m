function Y=nanmedfilt1(X,N,BlockSize)
% One dimensional median filter that ignores NaNs.
% Package: Util.cell
% Description: One dimensional median filter that ignores NaNs.
% Input  : - Vector on which to run the filter.
%          - Filter order, default is 3.
%            If order is odd then Y(k) is the std
%            of X( k-(N-1)/2 : k+(N-1)/2 ), else
%            it is the StD of X( k-N/2 : k+N/2-1 )
%          - BlockSize, by default its length of the input vector.
%            Use in case not enough computer memory.
% Output : - Median vector of the same length of the input vector.
%            Assuming the edges are equal to zero.
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Jun 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

Def.N = 3;
Def.BlockSize = length(X);

if (nargin==1)
   N         = Def.N;
   BlockSize = Def.BlockSize;
elseif (nargin==2)
   BlockSize = Def.BlockSize;
else
   error('Illegal number of input arguments');
end

NX = length(X);
if rem(N,2)~=1    % n even
    M = N/2;
else
    M = (N-1)/2;
end
Xe = [zeros(M,1); X; zeros(M,1)];
Y  = zeros(NX,1);

% Work in chunks to save memory
Indr = (0:N-1)';
Indc = 1:NX;
for I=1:BlockSize:NX
    Ind = Indc(ones(1,N),I:min(I+BlockSize-1,NX)) + ...
          Indr(:,ones(1,min(I+BlockSize-1,NX)-I+1));
    XX = reshape(Xe(Ind),N,min(I+BlockSize-1,NX)-I+1);
    Y(I:min(I+BlockSize-1,NX)) = nanmedian(XX,1);
end

