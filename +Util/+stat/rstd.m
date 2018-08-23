function Rstd=rstd(Mat,Dim)
% Robust std calculated from the 50% inner percentile of the data.
% Package: Util.stat
% Description: Robust std calculated from the 50% inner percentile
%              of the data.
% Input  : - Matrix.
%          - Dimension along to calculate the std. Default is 1.
% Output : - Robust std.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Rstd=rstd(randn(1000,3))
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1)
    Dim = 1;
end

Factor = 1.4826;  % = 1./norminv(0.75,0,1)

ValLow  = prctile(Mat,25,Dim);
ValHigh = prctile(Mat,75,Dim);
Rstd    = (ValHigh - ValLow).*0.5.*Factor;