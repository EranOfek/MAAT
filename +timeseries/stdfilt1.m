function Y=stdfilt1(X,N,BlockSize,Q)
% Standart deviation filter 
% Package: timeseries
% Description: One dimensional StD filter.
% Input  : - Vector on which to run the filter.
%          - Filter order, default is 3.
%            If order is odd then Y(k) is the std
%            of X( k-(N-1)/2 : k+(N-1)/2 ), else
%            it is the StD of X( k-N/2 : k+N/2-1 )
%          - BlockSize, by default its length of the input vector.
%            Use in case not enough computer memory.
%            If empty matrix then use default.
%          - If this parameter is provided then use quantile instaed
%            of std, where the parameter specify the fraction of data
%            to be in the returned range. For example 0.6834 is analog
%            to one sigma.
%            If two elements vector is given than these are the lower
%            and upper quantile range.
%            Default is empty matrix (i.e., use std).
% Output : - StD vector of the same length of the input vector.
%            Assuming the edges are equal to zero.
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Jun 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Y=timeseries.stdfilt1(rand(100,1),10,[],0.68);
% Reliable: 2
%--------------------------------------------------------------------------

Def.N = 3;
Def.BlockSize = length(X);
Def.Q         = [];
if (nargin==1)
   N         = Def.N;
   BlockSize = Def.BlockSize;
   Q         = Def.Q;
elseif (nargin==2)
   BlockSize = Def.BlockSize;
   Q         = Def.Q;
elseif (nargin==3)
   Q         = Def.Q;
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(BlockSize))
   BlockSize = Def.BlockSize;
end

if (~isempty(Q))
   if (length(Q)==1)
      Q = [(1-Q)./2];
      Q = [Q, 1-Q];
   end
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
    if (isempty(Q))
       % use std
       Y(I:min(I+BlockSize-1,NX)) = std(XX,0,1);
    else
       % use quantile
       Y(I:min(I+BlockSize-1,NX)) = quantile(XX,Q(2),1)-quantile(XX,Q(1),1);
    end
end

