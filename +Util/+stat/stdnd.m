function [Std]=stdnd(Data,Flag)
% Return the global StD of a N-D matrix.
% Package: Util.stat
% Description: Return the global StD of a N-D matrix.
% Input  : - N-D matrix.
%          - Flag: 0 to normalize by N-1; 1 to normalized by N,
%            default is 0.
% Output : - Value of global StD.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                      July 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: meannd.m, minnd.m, maxnd.m, mediannd.m
% Example: S=stdnd(rand(10,10,10),1)
% Reliable: 1
%------------------------------------------------------------------------------
if (nargin==1)
   Flag = 0;
end

%SizeData = size(Data);
%Vec = reshape(Data, [prod(SizeData), 1]);
%[Std] = std(Vec,Flag);

Std = std(Data(:),Flag);
