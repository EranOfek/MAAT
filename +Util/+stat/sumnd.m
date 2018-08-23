function Sum=sumnd(Data);
%------------------------------------------------------------------------------
% sumnd function                                                     AstroStat
% Description: Return the global sum of a N-D matrix.
% Input  : - N-D matrix
% Output : - Value of global sum.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                     March 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: meannd.m, minnd.m, maxnd.m, mediannd.m
% Example: Sum=sumnd(rand(4,4,4,4));
% Reliable: 1
%------------------------------------------------------------------------------

Sum = sum(Data(:));

%SizeData = size(Data);
%Vec = reshape(Data, [prod(SizeData), 1]);
%Sum = sum(Vec);


   
