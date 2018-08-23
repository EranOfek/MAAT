function Prct=prctile1(X,P)
% Faster version of prctile limited to vector input.
% Package: Util.math
% Description: calculate the percenile of a sample in a vector.
%              This is similar to the builtin function prctile.m but
%              limited to vector inputs and somewhat faster.
% Input  : - A vector.
%          - Percenile (in percents, i.e., 50).
% Output : - Percenile value.
% Tested : Matlab R2011b
%     By : Eran O. Ofek              Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Pr=Util.math.prctile1(rand(100,1),50);
% Reliable: 2
%---------------------------------------------------------------------------

P = P./100;
[N,M] = size(X);
if (min(N,M)>1)
   error('First argument for prcetile1 must be a vector');
end
N     = max(N,M);

SX   = sort(X);
Imin = floor(P.*N);
Imax = ceil(P.*N);
if (Imin<1)
   Imin = 1;
end
if (Imax<1)
   Imax = 1;
end
if (Imax==N)
   Imax = N;
end
if (Imin==N)
   Imin = N;
end
Prct = 0.5.*(SX(Imin) + SX(Imax));





