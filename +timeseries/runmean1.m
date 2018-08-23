function RM=runmean1(X,HalfSize)
% Running mean on equally spaced 1-D vector
% Package: timeseries
% Description: running mean on a 1-D vector.
% Input  : - Vector.
%          - Half size of running mean filter. Filter size will be
%            2*HalfSize+1. Default is 1.
% Output : - Output vector after running mean. The output vector will
%            have the same size as the input vector.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: R=zeros(10,1); R(5)=1; RM=timeseries.runmean1(R,2);
% Reliable: 1
%--------------------------------------------------------------------------


Def.HalfSize = 1;
if (nargin==1),
   HalfSize   = Def.HalfSize;
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end

B  = ones(2.*HalfSize+1,1)./(2.*HalfSize+1);
RM = conv(X,B,'same');
