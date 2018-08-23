function [ErrorMean,StD]=mean_error(M,Dim)
% Calculate the error on the mean using std/sqrt(N).
% Package: Util.stat
% Description: Calculate the error on the mean using std/sqrt(N).
% Input  : - An array.
%          - Diemsnion along to calculate the error on the mean.
%            Default is 1.
% Output : - Error on the mean.
%          - StD.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [ErrorMean,StD]=mean_error(rand(5,3,3),1)
% Reliable: 2
%--------------------------------------------------------------------------
StdFlag = 0;

if (nargin==1),
    Dim = 1;
end

StD       = squeeze(nanstd(M,StdFlag));
Nnn       = squeeze(sum(~isnan(M),Dim));
ErrorMean = StD./sqrt(Nnn);


