function [Median]=mediannd(Data);
%------------------------------------------------------------------------------
% mediannd function                                                  AstroStat
% Description: Return the global median of a N-D matrix.
% Input  : - N-D matrix
% Output : - Value of global Median.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                      July 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: meannd.m, minnd.m, maxnd.m, mediannd.m
% Example: M=mediannd(rand(10,10,10));
% Reliable: 1
%------------------------------------------------------------------------------

%SizeData = size(Data);
%Vec = reshape(Data, [prod(SizeData), 1]);
%[Median] = median(Vec);

Median = median(Data(:));


