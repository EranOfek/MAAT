function [rk,Q,S]=ccf_old(x1,x2,J)
% Cross coorelation between two equally spaced series
% Package: timeseries
% Description: Cross correlation function for two sets of equally
%              spaced time series. 
% Input  : - two column matrix in which the first column, is
%            an arbitrary index, and the second column is
%            the variable to correlate.
%          - second matrix of observations.
%	   - lag to run on (default is N/4).
% Output : - Two column matrix, in which the first column contains
%            the lag and the second the correlation.
%	   - Pormanteau statistics.
%	   - The standard deviation of the ccf (for very large N).
% Reference : Koen, C. & Lombard, F. 1993 MNRAS 263, 287-308
% See Also  : acf.m, ccf.m
% Tested : Matlab 4.2
%     By : Eran O. Ofek                    Nov 1995
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: timeseries.ccf_old([(1:1:100)',rand(100,1)],[(1:1:100)',rand(100,1)],20)
% Reliable: 2
%------------------------------------------------------------------------------
N1 = length(x1(:,1));
N2 = length(x2(:,1));
if N1~=N2,
   error('series must have the same length');
else
   N = N1;
end
c_x = 1;
c_y = 2;
if nargin==2,
   J   = floor(N/4);
elseif nargin==3,
   c_x = 1;
   c_y = 2;
   c_x = c_x;
else
   error('2 or 3 args only');
end
x1 = sortrows(x1,c_x);
x2 = sortrows(x2,c_x);
y1 = x1(:,c_y);
y2 = x2(:,c_y);
if J>N/2,
   'The lag is exceding size/2 of series'
end
my1  = mean(y1);
my2  = mean(y2);
MSy1 = y1 - my1;
MSy2 = y2 - my2;
Smy1 = std(MSy1);
Smy2 = std(MSy2);

for k=0:J,
   s = 0;
   for t=1:N-k,
      s = s + MSy1(t).*MSy2(t+k);
   end
   c(k+1) = (s./(N-k))./(Smy1.*Smy2);
end
Space = x1(2,c_x) - x1(1,c_x); 
rk = [Space.*rot90(0:J,3), rot90(c,3)];
Q = N.*sum(rk(1:J,2).*rk(1:J,2));
S = (1-rk(:,2).^2)./sqrt(N-1);
