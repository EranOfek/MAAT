function [NearestX,IndX,DeltaX]=nearest_unflag(X0,X,Flag)
% Nearest coordinate to an unflagged point.
% Package Util.array
% Description: Given a coordinate X0, a list of X-coordinates and
%              a flag (0 | 1) for each coordinate, find the nearest X
%              coordinate in the list to X0 that is Flag==0.
% Input  : - X0 coordinate for which to search the nearest coordinate.
%          - List of X coordinates.
%          - List of flags (0 | 1) which corresponds to the list of
%            coordinates, default is [0].
% Output : - Nearest coordinate tabulated in the list X for which
%            Flag==0 and its nearest to X0.
%          - Index of nearest coordinate.
%          - Absolute value difference between nearest coordinate and X0.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Mar 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [NearestX,Ind]=Util.array.nearest_unflag(0.74,rand(100,1),floor(rand(100,1).*2));
%-----------------------------------------------------------------------
if (nargin==2)
   Flag  = zeros(lensgth(X),1);
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end

I1 = find(Flag==0);
[DeltaX,MinInd] = min(abs((X(I1)-X0)));
IndX = I1(MinInd);
NearestX = X(IndX);
