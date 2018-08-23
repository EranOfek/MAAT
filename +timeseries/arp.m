function [a,cer,th]=arp(x,p,ex)
% Fit autoregression model for evently spaced time series
% Package: timeseries
% Description: Calculate the autoregessive process of order p.
%              Moddeling an evenly spaced time series, z(t), with:
%              z(t) = sum_i[a(i)*z(t-i)] + err(t) 
% Input  : - Matrix in which the first column is index (time)
%            and the second column is the observation value.
%	   - The autoregressive order (lag to run on). default is N/4.
%	   - The number of points beyond the last point of the series
%            to extraplate the series. Default is zero.
% Output : - AR(p) parameters.
%	   - AR(p) error in the parameters.
%	   - Series extrapolated into the future.
%	     The first column is index and the second is for the
%            predicted value of the series. 
% Reference : Koen, C. & Lombard, F. 1993 MNRAS 263, 287-308
% Tested : Matlab 4.2
%     By : Eran O. Ofek                    Sep 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [A,Cer]=timeseries.arp([(1:1:100)',rand(100,1)],10);
% Reliable: 2
%--------------------------------------------------------------------------

c_x = 1;
c_y = 2;
N = length(x(:,c_x));
if nargin==1,
   p   = floor(N/4);
elseif nargin==2,
   c_x = 1;
elseif nargin>3,
   error('1, 2 or 4 args only');
end
% x = sortby(x,c_x);
y = x(:,c_y);
if p>N/4,
   'The lag is exceding size/4 of series'
end
my = mean(y);
zt = y - my;
z = zeros(N-p,p);
k = 1;
for i=p+1:N,
   z(k,:) = rot90(zt(i-p:i-1),3);
   k = k + 1;
end
zts = zt(p+1:N,1);
a=z\zts;
cer = sqrt(diag(cov(z)));
if (nargin==3 & nargout==3),
   sp = x(2,c_x) - x(1,c_x);
   st = x(N,c_x) + sp;
   th = zeros(ex,2);
   for j=1:ex,
      th(j,c_x) = st + sp.*(j - 1);
      th(j,c_y) = rot90(a,1)*rot90(zt(i+j-p:i+j-1),2);
      zt(N+j)   = th(j,c_y);
   end
elseif (nargin==3 && nargout~=3),
   error('theoretical extension needed two output args')
end
